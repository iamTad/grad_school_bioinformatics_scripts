'''
This script is here to convert gff to fasta sequence
'''
import argparse
from Bio.SeqIO import FastaIO

def get_fasta(x):
    # takes in a fasta file and outputs a dictionary of sequences
    fa_dict = {}

    with open(x) as handle:
        for values in FastaIO.SimpleFastaParser(handle):
            fa_dict[values[0]] = values[1]

    return fa_dict


def get_gff(gff):
    # takes in gff file and outputs dictionary of {gene:{chrom: chromosome, strand: strand, range = [[start,end], [start,end]]}
    # assume gff file is correct/ordered
    gff_dict = {}

    with open(gff) as infile:
        for i in infile:
            if i.startswith('#'):
                continue
            line   = i.strip().split()
            chrom  = line[0]
            cat    = line[2]
            start  = int(line[3])
            end    = int(line[4])
            strand = line[6]
            gene   = line[-1].split(';')[0].split('-mRNA')[0][3:]

            if cat != 'CDS':
                continue

            if gene not in gff_dict:
                gff_dict[gene] = {'chrom': chrom,
                                  'strand': strand,
                                  'range': []}

            gff_dict[gene]['range'].append((start, end))

    return gff_dict


def main(argv):
    ### 'global' variables
    complement = {'A': 'T',
                  'T': 'A',
                  'C': 'G',
                  'G': 'C',
                  'N': 'N'}    

    ### initialize output
    ofile = open(argv.out, 'w')

    ### get fasta sequence
    fa_dict = get_fasta(argv.fas)

    ### get genes
    gene_dict = get_gff(argv.gff)

    ### start getting sequences
    for g in gene_dict:
        seq = ''
        c = gene_dict[g]['chrom']
        s = gene_dict[g]['strand']

        if argv.skip == True and c not in fa_dict:
            continue

        if s == '+':
            for r in gene_dict[g]['range']:
                seq += fa_dict[c][r[0]:r[1]+1]
        elif s == '-':
            for r in gene_dict[g]['range'][::-1]: # reverse order of ranges
                seq += ''.join([complement[i] for i in fa_dict[c][r[1]+1:r[0]:-1]])
        else:
            break
            print('strand is broken')
            return

        ofile.write(f'>{g}\n{seq}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fas', help = 'fasta file')
    parser.add_argument('-g', '--gff', help = 'gff file')
    parser.add_argument('-o', '--out', help = 'output file')
    parser.add_argument('-s', '--skip', help = 'skip nonmatching chromosomes (for sulfurigaster)', default = False, type = bool)
    args = parser.parse_args()

    main(args)
