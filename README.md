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

To make full use of the script you must set the environment variables:
`ENTREZ_EMAIL`
`ENTREZ_API_KEY`

You can generate your own API keys by creating an account on `https://www.ncbi.nlm.nih.gov/account/` with your ORCID.

```
pretextannotate -h

pretextannotate -i <input_file> -o <output_dir> -p <prefix> -c <context = '{"prim_accession"|"hap1_accession": "GCA_965178025.1"}'>
```

Context is a dictionary input containing the type of accession and the GCA accession number of the specific assembly.
