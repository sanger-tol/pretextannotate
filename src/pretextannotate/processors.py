import os
import json
import logging
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
from pretextannotate.chromosome_extraction import extract_chromosomes_only

logger = logging.getLogger('pretextannotation_logger')

def get_chromosomes(prefix: str, context: dict) -> list[dict] | None:
    """
    Get chromosome data for assembly
    """
    for key in ("prim_accession", "hap1_accession"):
        if key in context:
            return extract_chromosomes_only(context[key])

    logger.error(f"[Pretext Annotation] No accession in context for {prefix}")
    return None

def get_raw_length(context: dict, chroms: list[dict]) -> int | None:
    """
    Get raw length of genome
    """
    if context.get("genome_length_unrounded") or context.get("hap1_genome_length_unrounded"):
        return context.get("genome_length_unrounded") or context.get("hap1_genome_length_unrounded")
    else:
        logger.critical("""
    [Pretext Annotation] GENOME LENGTH NOT PROVIDED ON CLI
        - Keep in mind that the fall back method sums the lengths of chromosomes from the API
        - THIS SUM DOES NOT INCLUDE unlocs, if there are many then the plot may look off.
    """)
        return sum(chrom["length"] * 1e6 for chrom in chroms)

def fits_block(tw: float, block_width: float, fraction: float) -> bool:
    return tw <= block_width * fraction

def overlaps_prev(left: int, right: int, boxes, pad=0) -> bool:
    return any(left - pad < br and right + pad > bl for bl, br in boxes)

def add_mbp_scale(draw, font, left, top, w, h, total_length, font_size, text_colour):
    """Add Mbp scale to bottom of pretext map with smart positioning"""
    # Calculate reasonable tick intervals based on total genome size
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

    # Check for label crowding and adjust interval
    estimated_label_count: int = int(total_length / tick_interval) + 1
    avg_space_per_label: int = w / estimated_label_count if estimated_label_count > 0 else w
    sample_bbox = font.getbbox("1000")
    typical_label_width: int = sample_bbox[2] - sample_bbox[0]
    min_space_needed: int = typical_label_width * 1.5  # 50% padding between labels

    if avg_space_per_label < min_space_needed:
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

        logger.info(
            f"[Mbp Scale] Adjusted interval from original to {tick_interval} Mbp to prevent label crowding"
        )

    # Draw scale line
    scale_y: int = top + h + int(font_size * 0.5)
    draw.line([(left, scale_y), (left + w, scale_y)], fill=text_colour, width=2)

    # Pre-calculate the last tick that would get a label
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

def parse_context(context_dict_str: str) -> dict:
    try:
        return json.loads(context_dict_str)
    except Exception as e:
        logger.error(f"[Pretext Annotation] Failed to parse context_dict: {e}")
        return {}

def compute_chromosomes(prefix: str, context: dict, exclude: list[str], min_fraction: float) -> list[dict]:
    chroms = get_chromosomes(prefix, context)
    if chroms is None:
        raise SystemExit("Failed to retrieve chromosomes")
    max_length = max(c["length"] for c in chroms) if chroms else 0
    filtered = [
        c for c in chroms
        if c["molecule"] not in (exclude or [])
        and c["length"] >= min_fraction * max_length
    ]
    return sorted(filtered, key=lambda x: x["length"], reverse=True)

def choose_font_size(base_size: int, chrom_count: int, max_font_size: int = 90) -> int:
    if chrom_count > 25:
        return base_size
    if chrom_count > 20:
        return min(max_font_size, base_size + 5)
    if chrom_count > 10:
        return min(max_font_size, base_size + 10)
    return min(max_font_size, base_size + 20)

def build_canvas(image_path: str, font_path: str, font_size: int, background_colour: str):
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
    total = total_length or sum(c["length"] for c in sorted_chroms)

    acc, x_positions = 0, []
    for c in sorted_chroms:
        block = (c["length"] / total) * width
        x_positions.append(acc + block / 2)
        acc += block

    acc_h, y_positions = 0, []
    for c in sorted_chroms:
        block_h = (c["length"] / total) * height
        y_positions.append(acc_h + block_h / 2)
        acc_h += block_h

    return total, x_positions, y_positions

def draw_top_labels_with_positions(draw, font, sorted_chroms, total, left, font_size, text_colour, dot_width, max_fraction, x_positions, width):
    drawn_boxes = []
    for i, c in enumerate(sorted_chroms):
        label = str(c["molecule"])
        block: float = (c["length"] / total) * width
        bbox = font.getbbox(label)
        tw: float = bbox[2] - bbox[0]

        ok: bool = fits_block(tw, block, max_fraction)
        if ok:
            x_left: float = left + x_positions[i] - tw / 2
            x_right: float = x_left + tw
            if overlaps_prev(x_left, x_right, drawn_boxes):
                ok = False

        if ok:
            y = int(font_size * 0.6)
            draw.text((x_left, y), label, font=font, fill=text_colour)
            drawn_boxes.append((x_left, x_right))
        else:
            x: float = left + x_positions[i] - dot_width / 2
            y = int(font_size * 0.4)
            draw.text((x, y), ".", font=font, fill=text_colour)

def overlaps_prev_y(t, b, boxes, pad=0):
    return any(t - pad < bb and b + pad > tt for tt, bb in boxes)

def draw_left_labels(draw, font, sorted_chroms, total, left, top, height, font_size, text_colour, dot_width, vertical_label_field, y_positions, max_fraction):
    drawn_y_boxes: list[tuple[float, float]] = []

    for i, c in enumerate(sorted_chroms):
        lbl = str(c.get(vertical_label_field) or c.get("molecule") or "?")

        lb = font.getbbox(lbl)
        th: int = lb[3] - lb[1]
        tw: int = lb[2] - lb[0]

        block_h: int = int((c["length"] / total) * height)
        centre_y: int = top + y_positions[i]
        y_top: int = int(centre_y - th / 2 - font_size * 0.4)
        y_bot: int = y_top + th

        ok = block_h >= th * max_fraction
        if ok and overlaps_prev_y(y_top, y_bot, drawn_y_boxes, pad=5):
            ok = False

        if ok:
            x: float = left - tw - int(font_size * 0.4)
            draw.text((x, y_top), lbl, font=font, fill=text_colour)
            drawn_y_boxes.append((y_top, y_bot))
        else:
            dot_bbox = font.getbbox(".")
            dot_top: int = dot_bbox[1]
            dot_height: int = dot_bbox[3] - dot_bbox[1]
            x: float = left - dot_width - int(font_size * 0.4)
            y_dot: int = int(centre_y - (dot_height / 2) - dot_top)
            draw.text((x, y_dot), ".", font=font, fill=text_colour)

def label_pretext_map(args) -> tuple[Path, Path, Path]:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file_path: str = str(output_dir / f"{args.prefix}_annotated_pretext.png")

    context = parse_context(args.context_dict)

    logger.info(f"[Pretext Annotation] Input Snapshot Image is {args.pretext_file}")
    logger.info(f"[Pretext Annotation] Output will be saved at {output_file_path}")
    logger.info("[Pretext Annotation] Starting Pretext Annotation Process")

    sorted_chroms = compute_chromosomes(args.prefix, context, args.exclude_molecules, args.min_fraction)
    chrom_count = len(sorted_chroms)

    font_size = choose_font_size(args.font_size, chrom_count)
    logger.info(f"[Pretext Annotation] Adjusted font size: {font_size} for {chrom_count} chromosomes")

    canvas, draw, font, width, height, left, top = build_canvas(args.pretext_file, args.font, font_size, args.background_colour)

    dot_width: int = font.getbbox(".")[2] - font.getbbox(".")[0]
    raw_length = get_raw_length(context, sorted_chroms)
    total_length = (raw_length / 1e6) if raw_length else None

    logger.debug(f"[Pretext Annotation] total_length={total_length} Mb; chromosomes={chrom_count}")

    total, x_positions, y_positions = compute_positions(sorted_chroms, width, height, total_length)

    draw_top_labels_with_positions(draw, font, sorted_chroms, total, left, font_size, args.text_colour, dot_width, args.max_fraction, x_positions, width)
    draw_left_labels(draw, font, sorted_chroms, total, left, top, height, font_size, args.text_colour, dot_width, args.vertical_label_field, y_positions, args.max_fraction)

    if total_length:
        add_mbp_scale(draw, font, left, top, width, height, total_length, font_size, args.text_colour)

    canvas.save(output_file_path)
    logger.info(f"[Pretext Annotation] Saved labelled PNG → {output_file_path}")

    tif, gif = convert_png_to_tif_and_gif(str(output_file_path), dpi=(300, 300), max_width=1200)
    logger.info(f"[Pretext Annotation] Converted to TIFF & GIF → {tif}, {gif}")

    return Path(output_file_path), Path(tif), Path(gif)
