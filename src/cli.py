"""ToxProt CLI - Unified command-line interface for ToxProt data processing."""

import click

from src.data_processing.cli import data


@click.group()
@click.version_option(version="0.1.0", prog_name="toxprot")
def cli():
    """ToxProt - Animal toxin protein data processing toolkit.

    Process UniProt Swiss-Prot releases to extract toxin-related protein data.
    """
    pass


# Register command groups
cli.add_command(data)


if __name__ == "__main__":
    cli()
