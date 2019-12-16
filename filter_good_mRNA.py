'''
This script will check mRNA sequences for whether a transcript is a viable protein (from maker).
If it doesn't start with a start codon or has a premature stop codon or has no stop codon, then we
remove it from the fasta file.

Returns a fasta file
'''
from Bio.SeqIO import FastaIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help = 'input fasta file')
parser.add_argument('-o', '--out', help = 'output with only good sequences')
argv = parser.parse_args()

### Global Variables ###
codonStart = 'ATG'
codonStop  = set(['TAA', 'TAG', 'TGA'])

def main():
    outFile = open(argv.out, 'w')

    with open(argv.fasta) as handle:
        for values in FastaIO.SimpleFastaParser(handle):
            seqName  = values[0]
            sequence = values[1]
            ### check the sequences for "properness"
            # does it begin with a start codon?
            if not sequence.startswith(codonStart):
                continue
            # does it have a premature stop codon?
            for j in range(0,len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in codonStop and j != len(sequence) - 3:
                    # codon is a stop codon but it's not at the end, so it's premature
                    break
                ### updated part ###
                elif codon not in codonStop and j == len(sequence) -3:
                    break
                    
            # if we got this far, the sequence should be fine
            outFile.write('>{}\n{}\n'.format(seqName, sequence))

    return

if __name__ == '__main__':
    main()
