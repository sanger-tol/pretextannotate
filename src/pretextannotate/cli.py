import argparse
import pathlib
from pretextannotate.processors import label_pretext_map


def parse_args():
    parser = argparse.ArgumentParser(description="PretextAnnotate")

    # File Arguments
    parser.add_argument("--pretext-file", help="Input pretext PNG file", type=pathlib.Path)
    parser.add_argument("--prefix", help="Prefix for the output file", default="default")
    parser.add_argument("--context-dict", help="Context")
    parser.add_argument("--output", help="Output PNG file", default="./")

    # Font Arguments
    parser.add_argument("--font", help="Font file", default="./src/fonts/arial.ttf", type=pathlib.Path)
    parser.add_argument("--font-size", help="Font size", default=60)

    # Plot Arguments
    parser.add_argument("--exclude-molecules", help="List of molecules to exclude", default=[], type=list)
    parser.add_argument("--background-colour", help="Background colour", default="white")
    parser.add_argument("--text-colour", help="Text colour", default="black")
    parser.add_argument("--vertical-label-field", help="Vertical label field in output PNG", default="INSDC")

    # Other Arguments
    parser.add_argument("--min_fraction", help="Minimum Fraction of scaffolds to include", default="0.01")
    parser.add_argument("--max_fraction", help="Maximum Fraction of scaffolds to include", default="0.97")

    return parser.parse_args()

def main():
    args = parse_args()
    label_pretext_map(args)
