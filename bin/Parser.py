import csv

from Bio import SeqIO



def fasta_to_dict(fasta: str) -> dict:
    """read fasta and create dict with id and len of seq

    Args:
        fasta (str): path to fasta file

    Returns:
        dict: {'seq_id': len(seq)}
    """
    with open(fasta, "r") as f:
        seqs = {}
        for record in SeqIO.parse(f, "fasta"):
            seqs[record.id] = int(len(record.seq))
    return seqs


def slice_sequences(fasta: str, indexes: dict) -> dict:
    """slice each sequence using indexes

    Args:
        fasta (str): path to sequences to slice
        indexes (dict): dictionary of ids and indexes to slice

    Returns:
        dict: {id:seq}
    """
    seqs = {}
    with open(fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            id = record.id.split(":")[0]
            begin, end = indexes[id]
            # if end index is not in seq, take len(seq)
            if end > len(str(record.seq)):
                end = len(seq)
            seq = str(record.seq)[begin:end]
            # inf sliced string is not len(0), write to dict
            if not len(seq) == 0:
                seqs[id] = seq
    return seqs


def write_flanking(sliced_seqs: dict, outfile: str):
    """write fasta with new sequences

    Args:
        sliced_seqs (dict): {ids:(begin, end)}
        outfile (str): path to output file
    """
    with open("out.fasta", "w") as o:
        for k, v in sliced_seqs.items():
            o.write(f">{k.split('.')[0]}_flanking\n")
            o.write(f"{v}\n")


def blast(blast_result: str, fasta_dict: dict, flank: int) -> dict:
    """parse blast results and extract indexes of hits that meet the following requirements:

            qid = sid
            pid = 100.0
            length = q length
            mismatch = 0
            gaps = 0
        
        Flanking regions are included by taking start - 300 and end + 300. 

    Args:
        blast_result (str): path to blast results
        fasta_dict (dict): {id:len(seq)}

    Returns:
        dict: {id:(start, end)}
    """
    indexes = {}
    with open(blast_result) as br:
        reader = csv.reader(br, delimiter="\t")
        for line in reader:
            pass
            qid = line[0].split(":")[0]
            sid = line[1].split(":")[0]
            if (
                qid == sid
                and float(line[2]) == 100.0
                and int(line[3]) == fasta_dict[line[0]]
                and int(line[4]) == 0
                and int(line[5]) == 0
            ):
                sstart = int(line[8]) - int(flank)
                if sstart < 0:
                    sstart = 0
                ssend = int(line[9]) + int(flank)
                indexes[qid] = (sstart, ssend)
            else:
                print(f"{qid} did not satisfy requirements")

    return indexes
