#!/usr/bin/env python3
"""Convert SVG silhouettes to uniform grey PNGs."""

from io import BytesIO
from pathlib import Path

import cairosvg
from PIL import Image

# Configuration
SVG_DIR = Path("data/raw/PhyloPic/svg")
PNG_DIR = Path("data/raw/PhyloPic/png")

# Single uniform grey color for all silhouettes
GREY_VALUE = 102  # #666666


def svg_to_pil(svg_path: Path, scale: float = 2.0) -> Image.Image:
    """Convert SVG to PIL Image."""
    png_data = cairosvg.svg2png(url=str(svg_path), scale=scale)
    return Image.open(BytesIO(png_data)).convert("RGBA")


def convert_to_grey(img: Image.Image, alpha_threshold: int = 128) -> Image.Image:
    """Convert a silhouette image to uniform grey with binary alpha.

    Replaces all pixels above alpha threshold with fully opaque grey.
    Pixels below threshold become fully transparent.
    """
    img = img.copy()
    pixels = img.load()
    width, height = img.size

    for y in range(height):
        for x in range(width):
            *_, a = pixels[x, y]
            if a >= alpha_threshold:
                pixels[x, y] = (GREY_VALUE, GREY_VALUE, GREY_VALUE, 255)
            else:
                pixels[x, y] = (0, 0, 0, 0)

    return img


def main():
    """Convert SVG silhouettes to grey PNGs."""
    # Create output directory if needed
    PNG_DIR.mkdir(parents=True, exist_ok=True)

    # Special cases: use base files (Spider.svg, Scorpion.svg) instead of _black versions
    special_cases = {"Spider", "Scorpion"}

    # Find all *_black.svg files, excluding special cases
    svg_files = [
        f
        for f in sorted(SVG_DIR.glob("*_black.svg"))
        if f.stem.replace("_black", "") not in special_cases
    ]

    # Add special case files (without _black suffix)
    for name in special_cases:
        special_path = SVG_DIR / f"{name}.svg"
        if special_path.exists():
            svg_files.append(special_path)

    if not svg_files:
        print(f"No SVG files found in {SVG_DIR}")
        return

    print(f"Converting {len(svg_files)} SVG files to grey PNGs...")

    for svg_path in svg_files:
        # Generate output name: Cobra_black.svg -> Cobra_grey.png
        base_name = svg_path.stem.replace("_black", "")
        output_path = PNG_DIR / f"{base_name}_grey.png"

        print(f"  {svg_path.name} -> {output_path.name}")

        # Convert SVG to PIL Image, then to grey
        img = svg_to_pil(svg_path)
        grey_img = convert_to_grey(img)
        grey_img.save(output_path, "PNG")

    print(f"\nDone! Created {len(svg_files)} grey PNG files in {PNG_DIR}")


if __name__ == "__main__":
    main()
