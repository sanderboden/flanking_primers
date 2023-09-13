import argparse
import os

from Bio import Entrez

import Parser
import RunBlast



def cli():
    parser = argparse.ArgumentParser(
        prog='flanking-seq',
        description='extend DNA sequences with flanking regions',
        epilog='Created by Sander Boden (s.boden1@avans.nl)'
    )
    parser.add_argument(
        '-i', '--input', action='store',
        help='path to input fasta file'
    )
    parser.add_argument(
        '-o', '--output', action='store',
        help='path to output fasta',
        default='output.fasta'
    )
    parser.add_argument(
        '-l', '--length', action='store',
        help='length of flanking region to add',
        default=150
    )
    parser.add_argument(
        '-wd', '--working-dir', action='store',
        help='path to working directory to store intermediate results',
        default='workingdir'
    )
    parser.add_argument(
        '--api', action='store',
        help='NCBI api key',
        required=True
    )
    parser.add_argument(
        '--email', action='store',
        help='NCBI email-adress linked to api key',
        required=True
    )
    
    args = vars(parser.parse_args())
    
    
    if not os.path.exists(args['wd']):
        os.makedirs(args['wd'])
    else:
        print('working directory already exists. Exiting...')
        pass
    
    
    Entrez.api_key = args['api']
    Entrez.email = args['email']
    
    fasta = Parser.fasta_to_dict(args['i'])
    RunBlast.retrieve_contigs(fasta, f'{args["wd"]}/tmp.fasta')
    RunBlast.makedb(f'{args["wd"]}/tmp.fasta')
    RunBlast.run_blast(args['i'], f'{args["wd"]}/tmp.fasta')
    indexes = Parser.blast(f"{args['wd']}/result.txt", fasta)
    slices = Parser.slice_sequences(f'{args["wd"]}/tmp.fasta', indexes)
    Parser.write_flanking(slices, args['o'])



if __name__ == "__main__":
    cli()