import os
import json
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

def file_validator(in_file):
    """
    Validate files exist and are infact files
    """
    file_path = pathlib.Path(in_file)
    if file_path.exists() and file_path.is_file():
        logger.info(f"[Pretext Annotation] Input file valid: {file_path}")
        return file_path
    else:
        raise argparse.ArgumentTypeError(f"{in_file} is not a valid file")


def check_args(args):
    """
    More complex arg validation
    - context_dict and sizes are semi-mutually exclusive
        - if context_dict["genome_length"] and sizes exist, remove genome_length from the context_dict
        - if context_dict["accession"] does not exist then raise an error
        - if not context_dict or sizes exist, raise an error
    """
    context = json.loads(args.context_dict)

    if context and args.sizes:
        if args.sizes and context.get("genome_length"):
            data = context.get("genome_length")
            context.pop("genome_length")

            args.__dict__.update({"context_dict": f"{context}"})

            logger.error(f"[check_args] Can't provide genome length info ({data}) in --context AND provide a --sizes file! - This has been removed!")

    if context.get("accession"):
        logger.info(f"[check_args] Accession: {context['accession']}")
    else:
        raise argparse.ArgumentTypeError("[check_args] ACCESSION IS NEEDED")

    if args.context_dict is None and args.sizes is None:
        raise argparse.ArgumentTypeError("[check_args] Either context dict or sizes file is required")

    return args

def parse_args():
    parser = argparse.ArgumentParser(description="PretextAnnotate")

    # File Arguments
    parser.add_argument("-f", "--pretext_file", help="Input pretext PNG file", type=file_validator, required=True)
    parser.add_argument("-p", "--prefix", help="Prefix for the output file", default="default")
    parser.add_argument("-o", "--output", help="Output PNG file", default="./")

    ## Mutually exclusive
    parser.add_argument("-c", "--context_dict", help="Context")
    parser.add_argument("-s", "--sizes", help="Sizes file describing the input genome", type=file_validator)

    # Font Arguments
    data_path = os.path.join(os.path.dirname(__file__), 'fonts', 'OpenSans-Regular.ttf')
    logger.info(f"[Pretext Annotation] Found font file at {data_path}")

    parser.add_argument("--font", help="Font file", default=data_path, type=file_validator)
    parser.add_argument("--font_size", help="Font size", default=60, type=int)

    # Plot Arguments
    parser.add_argument("--exclude_molecules", help="List of molecules to exclude", nargs='+')
    parser.add_argument("--background_colour", help="Background colour", default="white")
    parser.add_argument("--text_colour", help="Text colour", default="black")
    parser.add_argument("--vertical_label_field", help="Vertical label field in output PNG", default="INSDC")

    # Other Arguments
    parser.add_argument("--min_fraction", help="Minimum Fraction of scaffolds to include", default=0.01, type=float)
    parser.add_argument("--max_fraction", help="Maximum Fraction of scaffolds to include", default=0.97, type=float)

    parser.add_argument("-v", "--versions", help="Return the version of the tool", action="version", version="%(prog)s: 1.1.0")

    return check_args(parser.parse_args())

def main():

    logger.info("[Pretext Annotation] Starting Pretext Annotation")

    args = parse_args()

    logger.info(f"[Pretext Annotation] PretextSnapshot: {args.pretext_file} | WITH | context_dict: {args.context_dict}")

    label_pretext_map(args)
