'''
This script will slide along the genome and give me the percent repeats as we
travel along the genome at some predetermined window size

output will be:
1. scaffold
2. midpoint
3. percent repeat
'''
##### modules #####
import argparse
from Bio.SeqIO import FastaIO


##### arguments #####
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help = 'fasta file input')
parser.add_argument('-w', '--window', help = 'window size to look for', default = 10000, type = int)
parser.add_argument('-m', '--missing', help = 'missing base value; is it "N"?', default = 'N')
parser.add_argument('-o', '--out', help = 'output file')
argv = parser.parse_args()


##### functions #####
def get_genome(x):
    '''
    input  = fasta file
    output = dictionary of fasta file
    '''
    fastaDict = {}

    with open(x) as handle:
        for values in FastaIO.SimpleFastaParser(handle):
            fastaDict[values[0]] = values[1]

    return fastaDict


def main():
    oFile = open(argv.out, 'w')

    fastaDict = get_genome(argv.fasta)

    for scaf in fastaDict:
        windowNumber = 1
        totalBins    = len(fastaDict[scaf])/argv.window # we're only using this to see if we continue the loop

        while windowNumber <= totalBins:
            sequence = fastaDict[scaf][(windowNumber-1)*argv.window:windowNumber*argv.window]
            #midpoint = (windowNumber-1)*argv.window+len(sequence)/2
            percRep  = sequence.count(argv.missing)/len(sequence)

            #oFile.write('{}\t{}\t{}\n'.format(scaf, int(midpoint), percRep))
            oFile.write('{}\t{}\t{}\t{}\n'.format(scaf, (windowNumber-1)*argv.window, (windowNumber-1)*argv.window + len(sequence), percRep))

            windowNumber += 1

    return 'done'

if __name__ == '__main__':
    main()
