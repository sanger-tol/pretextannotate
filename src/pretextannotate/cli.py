import os
import logging
import argparse
import pathlib
from pretextannotate.processors import label_pretext_map

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pretextannotation.log"),  # logs to file
        logging.StreamHandler()  # logs to console
    ]
)
logger = logging.getLogger('pretextannotation_logger')


def parse_args():
    parser = argparse.ArgumentParser(description="PretextAnnotate")

    # File Arguments
    parser.add_argument("--pretext_file", help="Input pretext PNG file", type=pathlib.Path)
    parser.add_argument("--prefix", help="Prefix for the output file", default="default")
    parser.add_argument("--context_dict", help="Context")
    parser.add_argument("--sizes", help="Sizes file describing the input genome", type=pathlib.Path)
    parser.add_argument("--output", help="Output PNG file", default="./")

    # Font Arguments
    data_path = os.path.join(os.path.dirname(__file__), 'fonts', 'OpenSans-Regular.ttf')
    logger.info(f"[Pretext Annotation] Found font file at {data_path}")
    parser.add_argument("--font", help="Font file", default=data_path, type=pathlib.Path)
    parser.add_argument("--font_size", help="Font size", default=60, type=int)

    # Plot Arguments
    parser.add_argument("--exclude_molecules", help="List of molecules to exclude", nargs='+')
    parser.add_argument("--background_colour", help="Background colour", default="white")
    parser.add_argument("--text_colour", help="Text colour", default="black")
    parser.add_argument("--vertical_label_field", help="Vertical label field in output PNG", default="INSDC")

    # Other Arguments
    parser.add_argument("--min_fraction", help="Minimum Fraction of scaffolds to include", default=0.01, type=float)
    parser.add_argument("--max_fraction", help="Maximum Fraction of scaffolds to include", default=0.97, type=float)

    parser.add_argument("-v", "--versions", help="Return the version of the tool", action="version", version="%(prog)s: 1.0.1")

    return parser.parse_args()

def main():

    logger.info("[Pretext Annotation] Starting Pretext Annotation")

    args = parse_args()

    logger.info(f"[Pretext Annotation] PretextSnapshot: {args.pretext_file} | WITH | context_dict: {args.context_dict}")

    label_pretext_map(args)
