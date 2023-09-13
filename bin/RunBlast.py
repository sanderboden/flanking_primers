import time
import subprocess

from Bio import Entrez



def retrieve_contigs(fasta_dict: dict, outfile) -> None:
    """download full contigs from ids in dict

    Args:
        fasta_dict (dict): dictionary with EntrezIDs as keys
    """
    for id in fasta_dict.keys():
        try:
            handle = Entrez.efetch(db="nucleotide", id=id.split(":")[0], retmode="xml")
            time.sleep(2)
        except:
            time.sleep(20)
            handle = Entrez.efetch(db="nucleotide", id=id.split(":")[0], retmode="xml")
        record = Entrez.read(handle)
        with open(outfile, "a") as out:
            out.write(f">{id}\n{record[0]['GBSeq_sequence']}\n")


def makedb(infile: str):
    """make blast db

    Args:
        infile (str): path to fasta
    """
    makeblastdb_cmd = [
        "makeblastdb",  # The makeblastdb program.
        "-in",
        infile,  # Replace with your input sequence file.
        "-dbtype",
        "nucl",  # Database type: 'nucl' for nucleotide sequences, 'prot' for protein sequences.
    ]

    try:
        subprocess.run(
            makeblastdb_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        print("BLAST database creation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database: {e}")
        print(e.stderr)


def run_blast(q: str, d: str, wd: str):
    """run blast

    Args:
        q (str): path to query fasta
        d (str): path to db fasta
    """
    blast_cmd = [
        "blastn",  # Replace with the appropriate BLAST program (e.g., blastp, blastx).
        "-query",
        q,  # Replace with your query sequence file.
        "-db",
        d,  # Replace with your BLAST database file.
        "-out",
        f"{wd}/result.txt",  # Specify the output file.
        "-evalue",
        "1e-5",  # Set the E-value threshold.
        "-outfmt",
        "6",  # Specify the output format (tabular).
    ]

    try:
        subprocess.run(
            blast_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        print("BLAST completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        print(e.stderr)
    time.sleep(2)
