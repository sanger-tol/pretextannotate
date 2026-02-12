# pretextannotate

PretextAnnotate is a script originally written by Karen van Niekerk (Sanger GRIT) in order to add chromosome annotations to a PNG image produced by PretextSnapshot.

`/src/fonts/arial.ttf` - Arial font file used for text rendering.
Font file is taken from https://github.com/kavin808/arial.ttf

`original_scripts` - Contains the original scripts used in this project, this is simply for archiving and reference purposes.

## Installation
```
git clone https://github.com/sanger-tol/pretextannotate.git

cd pretextannotate

pip install .
```

## Usage

To make full use of the script (and politeness to NCBI) you must set the environment variables:
`ENTREZ_EMAIL`
`ENTREZ_API_KEY`

You can generate your own API keys by creating an account on `https://www.ncbi.nlm.nih.gov/account/` with your ORCID.

```
pretextannotate -h

pretextannotate \\
    --pretext_file src/tests/ilDryDodo1.1_normal_FullMap.png \\
    --output ./ \\
    --prefix HELLO \\
    --context_dict '{"hap1_accession": "GCA_965178025.1"}'
```

Context is a dictionary input containing the type of accession and the GCA accession number of the specific assembly.

## Expected output

With the input pretext snapshot:
![ilDryDodo1 pretext map](./src/tests/ilDryDodo1.1_normal_FullMap.png)

As well as the arguments used in the Usage section.

The output should be a PNG, gif and tif file resembling:
![ilDryDodo1 annotatedpretext map](./src/tests/ilDryDodo1_Pretext.png)

There will also be a pretextannotation.log file containing, in this case:
```
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Starting Pretext Annotation
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] PretextSnapshot: src/tests/ilDryDodo1.1_normal_FullMap.png | WITH | context_dict: {"hap1_accession": "GCA_965178025.1"}
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Input Snapshot Image is src/tests/ilDryDodo1.1_normal_FullMap.png
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Output will be saved at .//HELLO_annotated_pretext.png
2026-02-12 12:55:46,100 [INFO] [Pretext Annotation] Starting Pretext Annotation Process
2026-02-12 12:55:46,470 [INFO] [Pretext Annotation] Adjusted font size: 60 for 31 chromosomes
2026-02-12 12:55:47,355 [INFO] [Pretext Annotation] Saved labelled PNG → .//HELLO_annotated_pretext.png
2026-02-12 12:55:47,783 [INFO] [Pretext Annotation] Converted .//HELLO_annotated_pretext.png → .//HELLO_annotated_pretext.tif, .//HELLO_annotated_pretext.gif
2026-02-12 12:55:47,783 [INFO] [Pretext Annotation] Converted to TIFF & GIF → .//HELLO_annotated_pretext.tif, .//HELLO_annotated_pretext.gif
```
