#!/usr/bin/env python3
"""Generate protein embeddings via Biocentral API.

Supports batching, deduplication, resume, and length-sorting.
Adapted from protspace-suite's biocentral embedding module.

Usage:
    python -m src.data_processing.embed_sequences input.fasta -o output.h5
    python -m src.data_processing.embed_sequences input.fasta --model prot_t5
"""

import argparse
import logging
import time
import warnings
from pathlib import Path

import h5py
import numpy as np
from biocentral_api import BiocentralAPI, CommonEmbedder, batched
from tqdm import tqdm

logger = logging.getLogger(__name__)

MODEL_SHORT_KEYS: dict[str, str] = {
    "prot_t5": "ProtT5",
    "prost_t5": "ProstT5",
    "esm2_8m": "ESM_8M",
    "esm2_650m": "ESM2_650M",
    "esm2_3b": "ESM2_3B",
}

DEFAULT_MODEL = "prot_t5"
DEFAULT_BATCH_SIZE = 1000


def resolve_embedder(name: str) -> str:
    """Resolve a short alias to the full embedder name."""
    if name in MODEL_SHORT_KEYS:
        return CommonEmbedder[MODEL_SHORT_KEYS[name]].value
    try:
        return CommonEmbedder[name].value
    except KeyError:
        pass
    # Assume it's already a full model name
    return name


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file into a dict of {entry_id: sequence}."""
    seqs: dict[str, str] = {}
    entry = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                entry = line[1:].split()[0]
                seqs[entry] = ""
            elif entry:
                seqs[entry] += line
    return seqs


def load_existing_ids(h5_path: Path) -> set[str]:
    """Return dataset keys already present in an HDF5 file."""
    if not h5_path.exists():
        return set()
    with h5py.File(h5_path, "r") as f:
        return set(f.keys())


def save_embeddings(h5_path: Path, embeddings: dict[str, np.ndarray]) -> None:
    """Append embeddings to an HDF5 file (one dataset per protein)."""
    with h5py.File(h5_path, "a") as f:
        for protein_id, emb in embeddings.items():
            if protein_id not in f:
                f.create_dataset(protein_id, data=emb.astype(np.float32))


def embed_fasta(
    fasta_path: Path,
    h5_path: Path,
    embedder: str = DEFAULT_MODEL,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Path:
    """Embed sequences from a FASTA file via Biocentral API.

    Supports deduplication, length-sorting, batching, and resume
    (skips IDs already present in h5_path).

    Returns the path to the completed HDF5 file.
    """
    embedder_name = resolve_embedder(embedder)
    sequences = read_fasta(fasta_path)
    logger.info("Read %d sequences from %s", len(sequences), fasta_path)

    # Resume: skip already-embedded sequences
    existing_ids = load_existing_ids(h5_path)
    if existing_ids:
        logger.info("Found %d existing embeddings in %s", len(existing_ids), h5_path)
    remaining = {k: v for k, v in sequences.items() if k not in existing_ids}
    logger.info(
        "Remaining: %d (skipped %d already embedded)",
        len(remaining),
        len(sequences) - len(remaining),
    )

    if not remaining:
        logger.info("All sequences already embedded in %s", h5_path)
        return h5_path

    # Deduplicate (API rejects batches with duplicate sequences)
    seq_to_ids: dict[str, list[str]] = {}
    for pid, seq in remaining.items():
        seq_to_ids.setdefault(seq, []).append(pid)

    n_duplicates = len(remaining) - len(seq_to_ids)
    if n_duplicates:
        logger.info(
            "%d duplicate sequences across %d proteins → %d unique to embed",
            n_duplicates,
            len(remaining),
            len(seq_to_ids),
        )

    # One representative ID per unique sequence
    unique_seqs = {ids[0]: seq for seq, ids in seq_to_ids.items()}

    # Sort by length to reduce padding waste
    unique_ids = sorted(unique_seqs.keys(), key=lambda pid: len(unique_seqs[pid]))

    # Connect to Biocentral
    logger.info("Connecting to Biocentral (%s)...", embedder_name)
    api = BiocentralAPI(fixed_server_url="https://biocentral.rostlab.org")
    api = api.wait_until_healthy(max_wait_seconds=30)
    logger.info("Server is healthy")

    # Batch and embed
    api_batches = list(batched(unique_ids, batch_size_limit=batch_size))
    total_embedded = 0
    failed_batches = 0

    pbar = tqdm(total=len(remaining), desc="Embedding", unit="seq")

    for batch_idx, batch_ids in enumerate(api_batches):
        batch_seqs = {pid: unique_seqs[pid] for pid in batch_ids}

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=".*longer than the recommended.*",
                    category=UserWarning,
                )
                result = api.embed(
                    embedder_name=embedder_name,
                    sequence_data=batch_seqs,
                    reduce=True,
                ).run()

            if result is not None:
                emb_dict = result.to_dict()
                if emb_dict:
                    # Expand embeddings to all IDs sharing the same sequence
                    expanded: dict[str, np.ndarray] = {}
                    for rep_id, emb in emb_dict.items():
                        seq = unique_seqs[rep_id]
                        for pid in seq_to_ids[seq]:
                            expanded[pid] = emb
                    save_embeddings(h5_path, expanded)
                    total_embedded += len(expanded)
            else:
                failed_batches += 1
                logger.error("Batch %d/%d: API returned None", batch_idx + 1, len(api_batches))

            batch_protein_count = sum(len(seq_to_ids[unique_seqs[pid]]) for pid in batch_ids)
            pbar.update(batch_protein_count)

        except Exception:
            failed_batches += 1
            logger.exception("Batch %d/%d failed", batch_idx + 1, len(api_batches))
            batch_protein_count = sum(len(seq_to_ids[unique_seqs[pid]]) for pid in batch_ids)
            pbar.update(batch_protein_count)

        if batch_idx < len(api_batches) - 1:
            time.sleep(1)

    pbar.close()

    print(f"\nDone. Embedded {total_embedded:,} / {len(remaining):,} sequences.")
    if failed_batches:
        print(f"Failed batches: {failed_batches} (rerun to retry)")
    print(f"Output: {h5_path}")

    return h5_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate protein embeddings via Biocentral API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Available models: {', '.join(sorted(MODEL_SHORT_KEYS))}",
    )
    parser.add_argument("fasta", type=Path, help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output H5 file (default: same name as FASTA with .h5 extension)",
    )
    parser.add_argument(
        "-m",
        "--model",
        default=DEFAULT_MODEL,
        help=f"Embedder model (default: {DEFAULT_MODEL})",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Sequences per API request (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )

    h5_path = args.output or args.fasta.with_suffix(".h5")
    embed_fasta(args.fasta, h5_path, embedder=args.model, batch_size=args.batch_size)


if __name__ == "__main__":
    main()
