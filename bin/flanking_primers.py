import argparse
import os

from Bio import Entrez

import Parser
import RunBlast





def cli():
    parser = argparse.ArgumentParser(
        prog='flanking-primers',
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
        '-w', '--workingdir', action='store',
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
    
    
    if not os.path.exists(args['workingdir']):
        os.makedirs(args['workingdir'])
    else:
        print('working directory already exists.')
    
    
    Entrez.api_key = args['api']
    Entrez.email = args['email']
    
    fasta = Parser.fasta_to_dict(args['input'])
    RunBlast.retrieve_contigs(fasta, f'{args["workingdir"]}/tmp.fasta')
    RunBlast.makedb(f'{args["workingdir"]}/tmp.fasta')
    RunBlast.run_blast(args['input'], f'{args["workingdir"]}/tmp.fasta', args['workingdir'])
    indexes = Parser.blast(f"{args['workingdir']}/result.txt", fasta, args['length'])
    slices = Parser.slice_sequences(f'{args["workingdir"]}/tmp.fasta', indexes)
    Parser.write_flanking(slices, args['output'])



if __name__ == "__main__":
    cli()