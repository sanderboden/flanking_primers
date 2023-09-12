import Parser
import RunBlast
from Bio import Entrez

qfasta = "in.fasta"
dbfasta = "tmp.fasta"

Entrez.api_key = "NCBI_API_KEY"
Entrez.email = "EMAIL"


if __name__ == "__main__":
    fasta = Parser.fasta_to_dict(qfasta)
    RunBlast.retrieve_contigs(fasta)
    RunBlast.makedb(dbfasta)
    RunBlast.run_blast(qfasta, dbfasta)
    indexes = Parser.blast("result.txt", fasta)
    slices = Parser.slice_sequences("tmp.fasta", indexes)
    Parser.write_flanking(slices)
