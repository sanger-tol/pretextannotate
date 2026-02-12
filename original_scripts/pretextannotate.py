#!/usr/bin/env python

"""
PretextAnnotate
-----------------
Originally written by Karen van Niekerk (kh18)

Updated by Damon-Lee B Pointon (dp24 / DLBPointon)
    - Converted script into a Python package
-----------------
Minimal, shareable utilities for labeling Pretext maps.
Extracted from scripts/process_images.py with non-Pretext functions removed.

Notes:
- Uses Pillow (imported as `PIL`): install with `pip install Pillow`.
- Requires `process_chromosome_data.py` in the same scripts folder.
- Assumes Pretext source PNGs live in `SERVER_DATA/pretext_images`.
"""

import logging
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from process_chromosome_data import extract_chromosomes_only

SERVER_DATA = "/Users/kh18/server_data"


def convert_png_to_tif_and_gif(png_path, dpi=(300, 300), max_width=None):
    """
    Given /…/Fig_N.png, writes:
      • /…/Fig_N.tif  at `dpi`
      • /…/Fig_N.gif  (optionally down‐scaled to max_width)
    Returns (tif_path, gif_path).
    These image conversions are needed for JATS XML submission, which requires TIFF for high-quality print and GIF for web display.
    """
    base, _ = os.path.splitext(png_path)
    img = Image.open(png_path)

    # optional resize to max_width
    if max_width and img.width > max_width:
        ratio = max_width / img.width
        new_size = (max_width, int(img.height * ratio))
        img = img.resize(new_size, Image.LANCZOS)

    # TIFF
    tif_path = f"{base}.tif"
    img.save(tif_path, format="TIFF", dpi=dpi)

    # GIF
    gif_path = f"{base}.gif"
    img.save(gif_path, format="GIF")

    logging.info(f"Converted {png_path} → {tif_path}, {gif_path}")
    return tif_path, gif_path


def resolve_open_sans_font(env_var="GENOMENOTES_FONT") -> str | None:
    """Locate an *upright* Open Sans font file.
    Preference order:
      1. $GENOMENOTES_FONT (explicit override)
      2. ~/Library/Fonts/OpenSans*.ttf (Homebrew / Font Book user install) **excluding Italic**
      3. vendored project_root/assets/fonts/OpenSans*.ttf (non-italic if available)
    Falls back to any OpenSans*.ttf if only Italic is found.
    Returns a string path or None if nothing found.
    """
    # 1. explicit override
    p = os.environ.get(env_var)
    if p and Path(p).is_file():
        return p

    # helper to pick first non-italic from list of paths
    def pick_upright(paths):
        upright = [x for x in paths if "italic" not in x.name.lower()]
        if upright:
            return str(sorted(upright)[0])
        return str(sorted(paths)[0]) if paths else None

    # 2. user fonts (Homebrew install ends up here)
    user_fonts = Path.home() / "Library" / "Fonts"
    hits = list(user_fonts.glob("OpenSans*.ttf"))
    chosen = pick_upright(hits)
    if chosen:
        return chosen

    # 3. vendored copy in repo assets
    script_root = Path(__file__).resolve().parent.parent
    pkg_dir = script_root / "assets" / "fonts"
    if pkg_dir.is_dir():
        hits = []
        for pat in ("OpenSans*.ttf", "open-sans*.ttf", "OpenSans-Regular.ttf"):
            hits.extend(pkg_dir.glob(pat))
        chosen = pick_upright(hits)
        if chosen:
            return chosen

    return None


def add_mbp_scale(draw, font, left, top, w, h, total_length, font_size, text_colour):
    """Add Mbp scale to bottom of pretext map with smart positioning"""
    # Calculate reasonable tick intervals based on total genome size
    if total_length <= 50:
        tick_interval = 10  # Every 10 Mbp
    elif total_length <= 200:
        tick_interval = 25  # Every 25 Mbp
    elif total_length <= 500:
        tick_interval = 50  # Every 50 Mbp
    elif total_length <= 1000:
        tick_interval = 100  # Every 100 Mbp
    elif total_length <= 2000:
        tick_interval = 200  # Every 200 Mbp
    else:
        tick_interval = 500  # Every 500 Mbp for very large genomes

    # Check for label crowding and adjust interval
    estimated_label_count = int(total_length / tick_interval) + 1
    avg_space_per_label = w / estimated_label_count if estimated_label_count > 0 else w
    sample_bbox = font.getbbox("1000")
    typical_label_width = sample_bbox[2] - sample_bbox[0]
    min_space_needed = typical_label_width * 1.5  # 50% padding between labels

    if avg_space_per_label < min_space_needed:
        if tick_interval == 10:
            tick_interval = 25
        elif tick_interval == 25:
            tick_interval = 50
        elif tick_interval == 50:
            tick_interval = 100
        elif tick_interval == 100:
            tick_interval = 200
        elif tick_interval == 200:
            tick_interval = 500
        elif tick_interval == 500:
            tick_interval = 1000

        logging.info(
            f"[Mbp Scale] Adjusted interval from original to {tick_interval} Mbp to prevent label crowding"
        )

    # Draw scale line
    scale_y = top + h + int(font_size * 0.5)
    draw.line([(left, scale_y), (left + w, scale_y)], fill=text_colour, width=2)

    # Pre-calculate the last tick that would get a label
    current_pos = 0
    tick_count = 0
    last_labeled_tick = None

    while current_pos <= total_length:
        if tick_count % (1 if tick_interval >= 50 else 2) == 0:
            x_pos = left + (current_pos / total_length) * w
            label = f"{int(current_pos)}"
            bbox = font.getbbox(label)
            label_width = bbox[2] - bbox[0]
            last_labeled_tick = (x_pos, label, label_width, current_pos)

        current_pos += tick_interval
        tick_count += 1

    # Calculate "Mbp" unit label dimensions
    unit_label = "Mbp"
    unit_bbox = font.getbbox(unit_label)
    unit_width = unit_bbox[2] - unit_bbox[0]

    # Determine optimal "Mbp" position
    min_gap = 15  # Minimum gap between labels
    right_margin_start = left + w  # End of the actual image
    ideal_mbp_x = right_margin_start + 15  # 15 pixels into the right margin

    skip_last_label = False
    final_mbp_x = ideal_mbp_x

    if last_labeled_tick:
        last_x, _last_label, last_width, _last_value = last_labeled_tick
        last_label_right = last_x + last_width / 2

        if ideal_mbp_x - min_gap < last_label_right:
            alt_mbp_x = last_label_right + min_gap

            if alt_mbp_x + unit_width <= left + w + 20:  # Allow slight overflow into right margin
                final_mbp_x = alt_mbp_x
            else:
                skip_last_label = True
                final_mbp_x = ideal_mbp_x

    # Draw all tick marks and labels
    current_pos = 0
    tick_count = 0

    while current_pos <= total_length:
        x_pos = left + (current_pos / total_length) * w

        # Draw tick mark
        tick_height = font_size // 4
        draw.line([(x_pos, scale_y), (x_pos, scale_y + tick_height)], fill=text_colour, width=2)

        # Add label (every other tick for readability if interval is small)
        if tick_count % (1 if tick_interval >= 50 else 2) == 0:
            label = f"{int(current_pos)}"
            bbox = font.getbbox(label)
            label_width = bbox[2] - bbox[0]

            # Skip the last tick label if we determined it should be skipped
            is_last_labeled_tick = (
                last_labeled_tick and current_pos == last_labeled_tick[3] and skip_last_label
            )

            if not is_last_labeled_tick:
                draw.text(
                    (x_pos - label_width / 2, scale_y + tick_height + font_size // 6),
                    label,
                    font=font,
                    fill=text_colour,
                )

        current_pos += tick_interval
        tick_count += 1

    # Add "Mbp" unit label at the calculated position
    draw.text(
        (final_mbp_x, scale_y + tick_height + font_size // 6),
        unit_label,
        font=font,
        fill=text_colour,
    )


def label_pretext_map(
    tolid: str,
    context: dict,
    output_dir: str,
    font_path: str | None = None,
    font_size: int = 60,
    exclude_molecules: list = None,
    min_fraction: float = 0.01,
    background_colour: str = "white",
    text_colour: str = "black",
    vertical_label_field: str = "INSDC",
) -> tuple[Path, Path, Path]:
    """
    1) Find server_data/pretext_images/{tolid}*.png
    2) Extract chrom-list + genome length from context
    3) Draw labels
    4) Save Fig_3_Pretext.png/.tif/.gif in output_dir
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) locate source folder by absolute path
    if font_path is None:
        font_path = resolve_open_sans_font()
        if font_path is None:
            font_path = "/System/Library/Fonts/Supplemental/Arial.ttf"  # fallback
    pretext_dir = Path(SERVER_DATA) / "pretext_images"
    logging.debug(f"[Pretext] Looking in {pretext_dir} for files named {tolid}*.png")

    matches = list(pretext_dir.glob(f"{tolid}*.png"))
    if not matches:
        logging.error(f"[Pretext] No .png files found for {tolid} in {pretext_dir}")
        return None

    src_png = matches[0]
    logging.info(f"[Pretext] Found source PNG: {src_png}")

    # 2) build chromosome_list & total_length
    if "prim_accession" in context:
        chroms = extract_chromosomes_only(context["prim_accession"])
    elif "hap1_accession" in context:
        chroms = extract_chromosomes_only(context["hap1_accession"])
    else:
        logging.error(f"[Pretext] No accession in context for {tolid}")
        return None

    raw_len = context.get("genome_length_unrounded") or context.get(
        "hap1_genome_length_unrounded"
    )
    total_length = (raw_len / 1e6) if raw_len else None
    logging.debug(f"[Pretext] total_length={total_length} Mb; chromosomes={len(chroms)}")

    # 3) draw labels
    img = Image.open(src_png)
    w, h = img.size
    try:
        font = ImageFont.truetype(font_path, font_size)
    except IOError:
        logging.warning(f"[Pretext] Could not load font {font_path}, using default")
        font = ImageFont.load_default()

    max_label_fraction = 0.97
    dot_width = font.getbbox(".")[2] - font.getbbox(".")[0]

    max_len = max(c["length"] for c in chroms)
    filtered = [
        c
        for c in chroms
        if c["molecule"] not in (exclude_molecules or [])
        and c["length"] >= min_fraction * max_len
    ]
    sorted_chroms = sorted(filtered, key=lambda x: x["length"], reverse=True)

    # Dynamically scale font size based on label density
    base_font_size = font_size
    chrom_count = len(sorted_chroms)
    max_font_size = 90

    if chrom_count > 25:
        font_size = base_font_size
    elif chrom_count > 20:
        font_size = min(max_font_size, base_font_size + 5)
    elif chrom_count > 10:
        font_size = min(max_font_size, base_font_size + 10)
    else:
        font_size = min(max_font_size, base_font_size + 20)

    logging.info(f"[Pretext] Adjusted font size: {font_size} for {chrom_count} chromosomes")

    top = int(font_size * 2.5)
    left = int(font_size * 7)
    bottom = int(font_size * 2.5)
    right = int(font_size * 2.5)

    canvas = Image.new("RGB", (w + left + right, h + top + bottom), background_colour)
    canvas.paste(img, (left, top))
    draw = ImageDraw.Draw(canvas)

    acc = 0
    x_positions = []
    total = total_length or sum(c["length"] for c in sorted_chroms)
    for c in sorted_chroms:
        block = (c["length"] / total) * w
        x_positions.append(acc + block / 2)
        acc += block

    def fits_block(tw, block_width, fraction=max_label_fraction):
        return tw <= block_width * fraction

    def overlaps_prev(l, r, boxes, pad=0):
        return any(l - pad < br and r + pad > bl for bl, br in boxes)

    # top labels
    drawn_boxes = []
    for i, c in enumerate(sorted_chroms):
        label = str(c["molecule"])
        block = (c["length"] / total) * w
        bbox = font.getbbox(label)
        tw = bbox[2] - bbox[0]

        ok = fits_block(tw, block)
        if ok:
            x_left = left + x_positions[i] - tw / 2
            x_right = x_left + tw
            if overlaps_prev(x_left, x_right, drawn_boxes):
                ok = False

        if ok:
            y = int(font_size * 0.6)
            draw.text((x_left, y), label, font=font, fill=text_colour)
            drawn_boxes.append((x_left, x_right))
        else:
            x = left + x_positions[i] - dot_width / 2
            y = int(font_size * 0.4)
            draw.text((x, y), ".", font=font, fill=text_colour)

    # side labels (vertical axis)
    y_positions = []
    acc_h = 0
    for c in sorted_chroms:
        block_h = (c["length"] / total) * h
        y_positions.append(acc_h + block_h / 2)
        acc_h += block_h

    drawn_y_boxes: list[tuple[float, float]] = []

    def overlaps_prev_y(t, b, boxes, pad=0):
        return any(t - pad < bb and b + pad > tt for tt, bb in boxes)

    for i, c in enumerate(sorted_chroms):
        lbl = str(c.get(vertical_label_field) or c.get("molecule") or "?")

        lb = font.getbbox(lbl)
        th = lb[3] - lb[1]
        tw = lb[2] - lb[0]

        block_h = (c["length"] / total) * h
        centre_y = top + y_positions[i]
        y_top = int(centre_y - th / 2 - font_size * 0.4)
        y_bot = y_top + th

        ok = block_h >= th * max_label_fraction
        if ok and overlaps_prev_y(y_top, y_bot, drawn_y_boxes, pad=5):
            ok = False

        if ok:
            x = left - tw - int(font_size * 0.4)
            draw.text((x, y_top), lbl, font=font, fill=text_colour)
            drawn_y_boxes.append((y_top, y_bot))
        else:
            dot_bbox = font.getbbox(".")
            dot_top = dot_bbox[1]
            dot_height = dot_bbox[3] - dot_bbox[1]
            x = left - dot_width - int(font_size * 0.4)
            y_dot = int(centre_y - (dot_height / 2) - dot_top)
            draw.text((x, y_dot), ".", font=font, fill=text_colour)

    if total_length:
        add_mbp_scale(draw, font, left, top, w, h, total_length, font_size, text_colour)

    # 4) save + convert
    out_png = output_dir / "Fig_3_Pretext.png"
    canvas.save(out_png)
    logging.info(f"[Pretext] Saved labelled PNG → {out_png}")

    tif, gif = convert_png_to_tif_and_gif(str(out_png), dpi=(300, 300), max_width=1200)
    logging.info(f"[Pretext] Converted to TIFF & GIF → {tif}, {gif}")

    return out_png, Path(tif), Path(gif)


if __name__ == "__main__":
    # Minimal standalone call to label one pretext map
    label_pretext_map(
        tolid="ilDryDodo1",
        context={"hap1_accession": "GCA_965178025.1"},
        output_dir="/Users/kh18/Documents/July_catchup/Drymonia_dodonaea",
    )
