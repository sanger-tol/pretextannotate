import os
import json
import logging
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from pretextannotate.chromosome_extraction import extract_chromosomes_only

logger = logging.getLogger('pretextannotation_logger')

def parse_sizes(sizes_file: Path) -> list[dict]:
    chrom_counter, chrom_data = int(), list()

    with open(sizes_file, 'r') as in_file:
        for line in in_file:
            molecular_name, molecule_length = line.split()
            chrom_counter += 1
            chrom_data.append({"INSDC": molecular_name, "length": int(molecule_length) / 1e6, "molecule": chrom_counter})

    return chrom_data


def get_raw_length(sizes, chroms: list[dict]) -> int:
    """
    Get raw length of genome
    """

    if sizes:
        logger.critical("[Pretext Annotation] Using sizes file for lengths")
    else:
        logger.critical("""
        [Pretext Annotation] GENOME LENGTH NOT PROVIDED ON CLI
            - Keep in mind that the fall back method sums the lengths of chromosomes from the API
            - THIS SUM DOES NOT INCLUDE unlocs, if there are many then the plot may look off.
        """)

    return sum(chrom["length"] for chrom in chroms) * 1e6

def fits_block(tw: float, block_width: float, fraction: float) -> bool:
    return tw <= block_width * fraction

def overlaps_prev(left: int, right: int, boxes, pad=0) -> bool:
    return any(left - pad < br and right + pad > bl for bl, br in boxes)

def calculate_tick_interval(total_length: int) -> int:
    """
    Calculate tick interval based on total genome length.
    """
    match True:
        case _ if total_length <= 50:
            tick_interval = 10  # Every 10 Mbp
        case _ if total_length <= 200:
            tick_interval = 25  # Every 25 Mbp
        case _ if total_length <= 500:
            tick_interval = 50  # Every 50 Mbp
        case _ if total_length <= 1000:
            tick_interval = 100  # Every 100 Mbp
        case _ if total_length <= 2000:
            tick_interval = 200  # Every 200 Mbp
        case _:
            tick_interval = 500  # Every 500 Mbp for very large genomes
    return tick_interval

def calculate_tick_interval_spacing(tick_interval: int, avg_space:int, min_space: int) -> int:
    """
    Recalculate the tick interval if spacing is too small between labels
    """
    if avg_space < min_space:
        match tick_interval:
            case 10:
                tick_interval = 25
            case 25:
                tick_interval = 50
            case 50:
                tick_interval = 100
            case 100:
                tick_interval = 200
            case 200:
                tick_interval = 500
            case 500:
                tick_interval = 1000
            case _:
                logger.error("[Add Mbp Scale] Tick interval larger than 500 MBP")
                raise ValueError("Tick interval larger than 500 MBP")

    return tick_interval

def calculate_last_labeled_tick(font, left, w, total_length, tick_interval):
    """
    Pre-calculate the last tick that would get a label
    """
    current_pos = 0
    tick_count = 0
    last_labeled_tick = None

    while current_pos <= total_length:
        if tick_count % (1 if tick_interval >= 50 else 2) == 0:
            x_pos: int = left + (current_pos / total_length) * w
            label = f"{int(current_pos)}"
            bbox = font.getbbox(label)
            label_width: int = bbox[2] - bbox[0]
            last_labeled_tick = (x_pos, label, label_width, current_pos)

        current_pos += tick_interval
        tick_count += 1

    return last_labeled_tick

def add_mbp_scale(draw, font, left, top, w, h, total_length, font_size, text_colour):
    """Add Mbp scale to bottom of pretext map with smart positioning"""

    first_tick_interval = calculate_tick_interval(total_length)
    logger.info(
        f"[Mbp Scale] Calculate a reasonable tick_interval based on genome size: {first_tick_interval}"
    )

    # Check for label crowding and adjust interval
    estimated_label_count: int = int(total_length / first_tick_interval) + 1
    avg_space_per_label: int = w / estimated_label_count if estimated_label_count > 0 else w
    sample_bbox = font.getbbox("1000")
    typical_label_width: int = sample_bbox[2] - sample_bbox[0]
    min_space_needed: int = typical_label_width * 1.5  # 50% padding between labels

    tick_interval = calculate_tick_interval_spacing(first_tick_interval, avg_space_per_label, min_space_needed)

    logger.info(
        f"[Mbp Scale] Adjusted interval from {first_tick_interval} to {tick_interval} Mbp to prevent label crowding"
    )

    # Draw scale line
    scale_y: int = top + h + int(font_size * 0.5)
    draw.line([(left, scale_y), (left + w, scale_y)], fill=text_colour, width=2)

    # Pre-calculate the last tick that would get a label
    last_labeled_tick = calculate_last_labeled_tick(font, left, w, total_length, tick_interval)

    # Calculate "Mbp" unit label dimensions
    unit_label = "Mbp"
    unit_bbox = font.getbbox(unit_label)
    unit_width: int = unit_bbox[2] - unit_bbox[0]

    # Determine optimal "Mbp" position
    min_gap = 15  # Minimum gap between labels
    right_margin_start: int = left + w  # End of the actual image
    ideal_mbp_x: int = right_margin_start + 15  # 15 pixels into the right margin

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
    tick_height = font_size // 4

    while current_pos <= total_length:
        x_pos = left + (current_pos / total_length) * w

        # Draw tick mark
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

def convert_png_to_tif_and_gif(png_path: str, dpi=(300, 300), max_width=None):
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
        img = img.resize(new_size, Image.Resampling.LANCZOS)

    # TIFF
    tif_path = f"{base}.tif"
    img.save(tif_path, format="TIFF", dpi=dpi)

    # GIF
    gif_path = f"{base}.gif"
    img.save(gif_path, format="GIF")

    logger.info(f"[Pretext Annotation] Converted {png_path} → {tif_path}, {gif_path}")
    return tif_path, gif_path

def compute_chromosomes(prefix: str, chroms: list[dict], exclude: list[str], min_fraction: float) -> list[dict]:
    """
    Figure out the chromosomes we want names in the image,
    this is based off of size of the molecule as well as whether the user wants to exclude certain molecule.
    """
    max_length = max(c["length"] for c in chroms) if chroms else 0
    filtered = [
        c for c in chroms
        if c["molecule"] not in (exclude or [])
        and c["length"] >= min_fraction * max_length
    ]
    logger.info(f"[compute chromosomes] {len(filtered)} chromosomes to be labelled, removed {len(chroms) - len(filtered)} molecules")

    return sorted(filtered, key=lambda x: x["length"], reverse=True)

def choose_font_size(base_size: int, chrom_count: int, max_font_size: int = 90) -> int:
    """
    Choose the font size based on the chromosome count
    """
    match True:
        case _ if chrom_count > 25:
            return base_size
        case _ if chrom_count > 20:
            return min(max_font_size, base_size + 5)
        case _ if chrom_count > 10:
            return min(max_font_size, base_size + 10)
        case _:
            return min(max_font_size, base_size + 20)

def build_canvas(image_path: str, font_path: str, font_size: int, background_colour: str):
    """
    Build the canvas to add the pretextmap to
    """
    image = Image.open(image_path)
    width, height = image.size
    font = ImageFont.truetype(font_path, font_size)

    top = int(font_size * 2.5)
    left = int(font_size * 7)
    bottom = int(font_size * 2.5)
    right = int(font_size * 2.5)

    canvas = Image.new("RGB", (width + left + right, height + top + bottom), background_colour)
    canvas.paste(image, (left, top))
    draw = ImageDraw.Draw(canvas)
    return canvas, draw, font, width, height, left, top

def compute_positions(sorted_chroms: list[dict], width: int, height: int, total_length: float | None):
    """
    Compute positions of the chromosome blocks on the pretext PNG
    This will be used to draw labels for the chromosome blocks
    """
    acc,acc_h,  x_positions, y_positions = 0, 0, [], []
    for c in sorted_chroms:
        block = (c["length"] / total_length) * width
        block_h = (c["length"] / total_length) * height
        x_positions.append(acc + block / 2)
        y_positions.append(acc_h + block_h / 2)
        acc += block
        acc_h += block_h

        logger.info(f"[compute positions] For {c['INSDC']} ({c['molecule']}) | X position = {x_positions[-1]},\tY position = {y_positions[-1]}")

    return x_positions, y_positions

def draw_top_labels_with_positions(draw, font, sorted_chroms, total, left, font_size, text_colour, dot_width, max_fraction, x_positions, width):
    """
    Draw the labels for the top of the image.
    Using the chromosome positions to set top labels (molecule number)
    """
    y_label = int(font_size * 0.6)
    y_dot = int(font_size * 0.4)
    drawn_boxes = []

    for i, c in enumerate(sorted_chroms):
        label = str(c["molecule"])
        block: float = (c["length"] / total) * width
        bbox = font.getbbox(label)
        text_width: int = bbox[2] - bbox[0]

        ok = fits_block(text_width, block, max_fraction)
        if ok:
            x_left: float = left + x_positions[i] - text_width / 2
            x_right: float = x_left + text_width
            if overlaps_prev(x_left, x_right, drawn_boxes):
                ok = False

        if ok:
            draw.text((x_left, y_label), label, font=font, fill=text_colour)
            drawn_boxes.append((x_left, x_right))
        else:
            x_dot = left + x_positions[i] - dot_width / 2
            draw.text((x_dot, y_dot), ".", font=font, fill=text_colour)

def overlaps_prev_y(top, bottom, boxes, padding=0):
    return any(top - padding < bb and bottom + padding > tt for tt, bb in boxes)

def draw_left_labels(draw, font, sorted_chroms, total, left, top, height, font_size, text_colour, dot_width, vertical_label_field, y_positions, max_fraction):
    """
    Draw the labels for the left of the image.
    Using the chromosome positions to set Y labels (molecule name)
    """
    drawn_y_boxes: list[tuple[float, float]] = []
    y_pad = int(font_size * 0.4)
    dot_bbox = font.getbbox(".")
    dot_top: int = dot_bbox[1]
    dot_height: int = dot_bbox[3] - dot_bbox[1]

    for i, c in enumerate(sorted_chroms):
        lbl = str(c.get(vertical_label_field) or c.get("molecule") or "?")

        lb = font.getbbox(lbl)
        text_height: int = lb[3] - lb[1]
        text_width: int = lb[2] - lb[0]

        block_h: float = (c["length"] / total) * height
        centre_y: int = top + y_positions[i]
        y_top = int(centre_y - text_height / 2 - y_pad)
        y_bot: int = y_top + text_height

        ok: bool = block_h >= text_height * max_fraction
        if ok and overlaps_prev_y(y_top, y_bot, drawn_y_boxes, padding=5):
            ok = False

        if ok:
            x: int = left - text_width - y_pad
            draw.text((x, y_top), lbl, font=font, fill=text_colour)
            drawn_y_boxes.append((y_top, y_bot))
        else:
            x_dot: int = left - dot_width - y_pad
            y_dot = int(centre_y - (dot_height / 2) - dot_top)
            draw.text((x_dot, y_dot), ".", font=font, fill=text_colour)

def label_pretext_map(args) -> tuple[Path, Path, Path]:
    """
    Main function for adding labels to a pretextmap
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path: str = str(output_dir / f"{args.prefix}_annotated_pretext.png")

    logger.info(f"[Pretext Annotation] Input Snapshot Image is {args.pretext_file}")
    logger.info(f"[Pretext Annotation] Output will be saved at {output_file_path}")
    logger.info("[Pretext Annotation] Starting Pretext Annotation Process")

    if args.sizes:
        logger.info(f"[Pretext Annotation] Sizes file provided: {args.sizes}")
        chroms = parse_sizes(args.sizes)
    else:
        logger.info("[Pretext Annotation] Sizes file not provided, falling back to NCBI api")
        chroms: list[dict[str, str]] = extract_chromosomes_only(args.accession)

    if not chroms:
        raise ValueError(f"NO CHROMOSOMES FOUND: {chroms}")

    sorted_chroms = compute_chromosomes(args.prefix, chroms, args.exclude_molecules, args.min_fraction)
    chrom_count = len(sorted_chroms)
    raw_length = get_raw_length(args.sizes, sorted_chroms)
    if raw_length == 0:
        raise ValueError(f"NO LENGTHS OF MOLECULE FOUND IN CHROM_LIST: {sorted_chroms}")
    total_length: float = ( raw_length / 1e6 )
    logger.debug(f"[Pretext Annotation] total_length={total_length} Mb; chromosomes={chrom_count}")


    font_size = choose_font_size(args.font_size, chrom_count)
    logger.info(f"[Pretext Annotation] Adjusted font size: {font_size} for {chrom_count} chromosomes")

    canvas, draw, font, width, height, left, top = build_canvas(args.pretext_file, args.font, font_size, args.background_colour)
    x_positions, y_positions = compute_positions(sorted_chroms, width, height, total_length)
    dot_width: int = font.getbbox(".")[2] - font.getbbox(".")[0]

    draw_top_labels_with_positions(draw, font, sorted_chroms, total_length, left, font_size, args.text_colour, dot_width, args.max_fraction, x_positions, width)
    draw_left_labels(draw, font, sorted_chroms, total_length, left, top, height, font_size, args.text_colour, dot_width, args.vertical_label_field, y_positions, args.max_fraction)
    add_mbp_scale(draw, font, left, top, width, height, total_length, font_size, args.text_colour)

    canvas.save(output_file_path)
    logger.info(f"[Pretext Annotation] Saved labelled PNG → {output_file_path}")

    tif, gif = convert_png_to_tif_and_gif(str(output_file_path), dpi=(300, 300), max_width=1200)
    logger.info(f"[Pretext Annotation] Converted to TIFF & GIF → {tif}, {gif}")

    return Path(output_file_path), Path(tif), Path(gif)
