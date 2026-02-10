"""ProtSpace analysis module for ToxProt CLI.

This module provides tools for generating and visualizing protein space
embeddings using ProtT5-XL-UniRef50 and UMAP dimensionality reduction.
"""

from .cli import protspace
from .config import (
    DEFAULT_MIN_DIST,
    DEFAULT_N_NEIGHBORS,
    TOP_N,
    VARIANT_CONFIGS,
)

__all__ = [
    "protspace",
    "TOP_N",
    "DEFAULT_N_NEIGHBORS",
    "DEFAULT_MIN_DIST",
    "VARIANT_CONFIGS",
]
