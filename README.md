# flanking_primers
Extract flanking regions for a fasta file from NCBI.


## What is this?
Use this if you have a file of NCBI sequences (such as gene sequences) (with headers like 'NKKQ01000016.1:....') and want to extend these sequences with flanking regions of x bp for primer design.
flanking-seq downloads the full contigs from which the sequence originates and uses BLAST to find and extend the indexes of the smaller sequences, including flanking regions.


## How to install
First, make sure python (3.10) and blast are in your path.
The easiest way (until I make the process easier) is to run the script using a virtual env
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

after running the script, deactivate the environment: `deactivate`


## Usage

```bash
python3 bin/flanking_primers.py -h
usage: flanking-primers [-h] [-i INPUT] [-o OUTPUT] [-l LENGTH] [-wd WORKING_DIR] --api API --email EMAIL

extend DNA sequences with flanking regions

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input fasta file
  -o OUTPUT, --output OUTPUT
                        path to output fasta
  -l LENGTH, --length LENGTH
                        length of flanking region to add
  -wd WORKING_DIR, --working-dir WORKING_DIR
                        path to working directory to store intermediate results
  --api API             NCBI api key
  --email EMAIL         NCBI email-adress linked to api key

Created by Sander Boden (s.boden1@avans.nl)
```
