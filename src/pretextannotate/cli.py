import argparse
import logging
import os
import pathlib

from pretextannotate.processors import label_pretext_map

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("pretextannotation.log"),  # logs to file
        logging.StreamHandler(),  # logs to console
    ],
)
logger = logging.getLogger("pretextannotation_logger")


def file_validator(in_file, file_type):
    """
    Validate files exist and are infact files
    """
    file_checks = {
        "font": {"format": [".tff"], "validated": False},
        "index": {"format": [".sizes", ".fai", ".genome"], "validated": False},
        "pretext_png": {"format": [".png"], "validated": False},
    }
    file_path = pathlib.Path(in_file)
    if file_path.exists() and file_path.is_file():
        if file_type == "font" and file_path.suffix in file_checks["font"]["format"]:
            file_checks["font"]["validated"] = True
            return file_path
        if file_type == "index" and file_path.suffix in file_checks["index"]["format"]:
            file_checks["index"]["validated"] = True
            return file_path
        if file_type == "pretext_png" and file_path.suffix in file_checks["pretext_png"]["format"]:
            file_checks["pretext_png"]["validated"] = True
            return file_path

        logger.info(f"[Pretext Annotation] Input file check results: {file_checks[file_type]}")
        return file_path
    else:
        raise argparse.ArgumentTypeError(f"{in_file} is not a valid file")


def check_args(args):
    """
    Expandable function to validate arguments and/or give some info
    """

    if args.index is None:
        logger.info(
            "[check_args] Without the sizes file, falling back to NCBI API using the Accession."
        )

    return args


def parse_args():
    parser = argparse.ArgumentParser(description="PretextAnnotate")

    # File Arguments
    parser.add_argument(
        "-f",
        "--pretext_file",
        help="Input pretext PNG file",
        type=lambda s: file_validator(s, "pretext_png"),
        required=True,
    )
    parser.add_argument("-p", "--prefix", help="Prefix for the output file", default="default")
    parser.add_argument("-o", "--output", help="Output PNG file", default="./")
    parser.add_argument(
        "-i",
        "--index",
        help="Index file describing the input genome",
        type=lambda s: file_validator(s, "index"),
    )

    # Font Arguments
    data_path = os.path.join(os.path.dirname(__file__), "fonts", "OpenSans-Regular.ttf")
    logger.info(f"[Pretext Annotation] Found font file at {data_path}")

    parser.add_argument(
        "--font",
        help="Font file",
        default=data_path,
        type=lambda s: file_validator(s, "font"),
    )
    parser.add_argument("--font_size", help="Font size", default=60, type=int)

    # Plot Arguments
    parser.add_argument("--exclude_molecules", help="List of molecules to exclude", nargs="+")
    parser.add_argument("--background_colour", help="Background colour", default="white")
    parser.add_argument("--text_colour", help="Text colour", default="black")
    parser.add_argument(
        "--vertical_label_field", help="Vertical label field in output PNG", default="INSDC"
    )

    # Other Arguments
    parser.add_argument(
        "--gca_accession", help="The GCA_Accession of the sample of interest - optional", type=str
    )
    parser.add_argument(
        "--min_fraction", help="Minimum Fraction of scaffolds to include", default=0.01, type=float
    )
    parser.add_argument(
        "--max_fraction", help="Maximum Fraction of scaffolds to include", default=0.97, type=float
    )

    parser.add_argument(
        "-v",
        "--version",
        help="Return the version of the tool",
        action="version",
        version="%(prog)s: 1.1.3",
    )

    return check_args(parser.parse_args())


def main():
    logger.info("[Pretext Annotation] Starting Pretext Annotation")

    args = parse_args()

    logger.info(f"[Pretext Annotation] PretextSnapshot: {args.pretext_file}")

    label_pretext_map(args)
