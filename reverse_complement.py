'''
This takes in a fasta file and a list of contigs to reverse complement.
Spits out a new fasta file that reverses the contigs of interest
'''
from Bio.SeqIO import FastaIO
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-fi', '--ifasta', help = 'input fasta file')
parser.add_argument('-r', '--rev', help = 'list of contigs to reverse')
parser.add_argument('-fo', '--ofasta', help = 'output fasta file')
argv   = parser.parse_args()


def rev_comp(x):
    # x is the sequence
    complementDict = {'C': 'G',
                      'G': 'C',
                      'A': 'T',
                      'T': 'A',
                      'N': 'N'}

    revComp = ''.join([complementDict[i.upper()] for i in x[::-1]])

    return revComp


def get_rc_list(x):
    # x is the reverse complement list
    rcList = []
    with open(x, 'r') as infile:
        for i in infile:
            line = i.strip()
            rcList.append(line)

    rcSet = set(rcList)

    return rcSet


def main():
    rcSet = get_rc_list(argv.rev)
    outfile = open(argv.ofasta, 'w')
    
    with open(argv.ifasta) as handle:
        for values in FastaIO.SimpleFastaParser(handle):
            contig  = values[0]
            if contig in rcSet:
                seq = rev_comp(values[1])
            else:
                seq = values[1]

            outfile.write('>{}\n{}\n'.format(contig, seq))

    return 'im pretty tired of this shit'


if __name__ == '__main__':
    main()
