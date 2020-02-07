'''
Write an HKA test for things, because I can't find one online.
Assume sites are biallelic

Filter out vcf for individuals such that at least 66% of individuals per population exists

aside from vcf input, we want population input:
- for each population/species, provide a file with one individual per line (col 1) and the population name
- each individual must be labeled as such in the vcf file
- population A will be the first unique occurring population
- population B will be the second unique occurring population
- make the outgroup the last unique appearing population

for HKA:
W = polymorphic in species A
X = polymorphic in species B
Y = fixed A and different from B and C
Z = fixed B and different from A and C 


Output should be:
1. chrom
2. start
3. end
4. hka_A
5. hka_B
6. homogeneity
'''
### Modules ###
import argparse
import io
import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


### Arguments ###
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', help = 'input vcf file')
parser.add_argument('-p', '--pop', help = 'list of populations')
parser.add_argument('-w', '--window', help = 'window size; assume step is same', type = int, default = 10000)
parser.add_argument('-c', '--cutoff', help = 'how many missing individuals max?', type = float, default = 0.4)
parser.add_argument('-o', '--out', help = 'output file name')
argv = parser.parse_args()


### Functions ###
def parse_populations(popFile):
    ### input is the population file
    ### output a dictionary of individual:species
    popLabelDict = {}
    popList      = []
    with open(popFile, 'r') as inFile:
        for i in inFile:
            line       = i.strip().split()
            sample     = line[0]
            population = line[1]

            if population not in popLabelDict:
                popLabelDict[population] = []
                popList.append(population)
            popLabelDict[population].append(sample)

    return popLabelDict, popList


def read_vcf(vcf):
    ### inputs a vcf file
    ### outputs a pandas dataframe
    with open(vcf, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]

    return pd.read_csv(
        io.StringIO(''.join(lines[1:])),
        names=lines[0].strip().split(),
        dtype={'#CHROM': str, 'POS': int},
        sep='\t').rename(columns={'#CHROM': 'CHROM'}), lines[0].strip().split()[9:]

def process_vcf_for_hka(df, head, pop_label):
    ### inputs a pandas dataframe
    ### outputs a dataframe with WXYZ
    
    vcfDF        = df
    header       = head
    popLabelDict = pop_label
    #print(header)
    #vcfDF.to_csv(r'/global/scratch/dmai/projects/phdPaper2/selectionScan/ank/hka/vcfdf.txt', index = False)
    # get genotype data for each sample
    for indiv in header:
        vcfDF[indiv] = vcfDF[indiv].str[0:1] + vcfDF[indiv].str[2:3]
    
    # concatenate samples of same species and get missing
    for i in popLabelDict:
        #vcfDF[i] = np.add.reduce(vcfDF[popLabelDict[i]], axis=1)
        vcfDF[i] = vcfDF[popLabelDict[i]].apply(lambda row: ''.join(row.values.astype(str)), axis=1)

    # remove extraneous columns including general data and samples
    vcfDF = vcfDF.drop(vcfDF.columns[range(2,9)], axis=1)
    vcfDF.drop(header, axis=1, inplace=True)
    
    # remove sites that are missing more than 70% individuals
    for i in popLabelDict:
        vcfDF = vcfDF.drop(vcfDF[vcfDF[i].str.count('\.')/vcfDF[i].str.len() > argv.cutoff].index)
    
    # make a set of the different alleles
    for i in popLabelDict:
        vcfDF[i] = vcfDF[i].apply(lambda x: ''.join(set([j for j in x if j != '.'])))
    
    return vcfDF

def main():
    outFile       = open(argv.out, 'w')
    popLabelDict, popList  = parse_populations(argv.pop)
    vcfDF, header = read_vcf(argv.vcf)
    vcfProcessed  = process_vcf_for_hka(vcfDF, header, popLabelDict)

    # unpack populations
    A, B, C = popList

    ### these are examples of how to further process my dataframe cuz I'm going to sleep
    # to filter by position: stuffC = stuffB[(stuffB['POS'] >= 10000) & (stuffB['POS'] <= 10200)]
    # to get a set of alleles: stuffB['A'] = stuffB['A'].apply(lambda x: ''.join(set([i for i in x if i != '.'])))
    # to subsample a dataframe: stuffB[stuffB['A'] == '1']['A_missing'].max()
    
    # Moving forward, I hardcode ABC, because I don't know how to do the inequalities for YZ in a generalized manner
    # get global values of WXYZ
    wGlobal = vcfProcessed[vcfProcessed[A].str.len() > 1].shape[0] # number of rows of polymorphic A
    xGlobal = vcfProcessed[vcfProcessed[B].str.len() > 1].shape[0] # number of rows of polymorphic B
    yGlobal = vcfProcessed[(vcfProcessed[A].str.len() == 1) &
                            (vcfProcessed[B].str.len() == 1) &
                            (vcfProcessed[C].str.len() == 1) &
                            (vcfProcessed[A] != vcfProcessed[B]) &
                            (vcfProcessed[B] == vcfProcessed[C])].shape[0] # number of rows of fixed A != B == C
    zGlobal = vcfProcessed[(vcfProcessed[A].str.len() == 1) &
                            (vcfProcessed[B].str.len() == 1) &
                            (vcfProcessed[C].str.len() == 1) &
                            (vcfProcessed[A] != vcfProcessed[B]) &
                            (vcfProcessed[A] == vcfProcessed[C])].shape[0] # number of rows of fixed B != A == C

    # loop through muller elements
    for c in vcfProcessed.CHROM.unique():
        # initialize window parameters
        windowCounter = 0
        # I just don't want to write a long expression in the while loop, so I make these two and will update in loop as well
        windowMin     = argv.window * windowCounter
        windowMax     = windowMin + argv.window        
        while windowMin <= vcfProcessed[vcfProcessed['CHROM'] == c]['POS'].max():
            # make a tiny vcf to work with so I don't keep specifying chromosome and range for local wxyz values
            miniDF = vcfProcessed[(vcfProcessed['CHROM'] == c) & 
                                  (vcfProcessed['POS'] >= windowMin) &
                                  (vcfProcessed['POS'] < windowMax)]

            # get local WXYZ values
            wLocal = miniDF[miniDF[A].str.len() > 1].shape[0]
            xLocal = miniDF[miniDF[B].str.len() > 1].shape[0]
            yLocal = miniDF[(miniDF[A].str.len() == 1) &
                            (miniDF[B].str.len() == 1) &
                            (miniDF[C].str.len() == 1) &
                            (miniDF[A] != miniDF[B]) &
                            (miniDF[B] == miniDF[C])].shape[0]
            zLocal = miniDF[(miniDF[A].str.len() == 1) &
                            (miniDF[B].str.len() == 1) &
                            (miniDF[C].str.len() == 1) &
                            (miniDF[A] != miniDF[B]) &
                            (miniDF[A] == miniDF[C])].shape[0]

            # get HKA for population A
            if wLocal + yLocal > 0:
                Ahka = -np.log10(chi2_contingency(np.array([[wLocal, yLocal], [wGlobal, yGlobal]]))[1])
            else:
                Ahka = np.nan

            # get HKA for population B
            if xLocal + zLocal > 0:
                Bhka = -np.log10(chi2_contingency(np.array([[xLocal, zLocal], [xGlobal, zGlobal]]))[1])
            else:
                Bhka = np.nan

            # get homogeneity score
            if (wLocal + yLocal > 0) and (wLocal + xLocal > 0) and (xLocal + zLocal > 0) and (yLocal + zLocal > 0):
                homo = -np.log10(chi2_contingency(np.array([[wLocal, yLocal], [xLocal, zLocal]]))[1])
            else:
                homo = np.nan

            outFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, windowMin, windowMax, 
                                                                            wLocal, xLocal, yLocal, zLocal,
                                                                            Ahka, Bhka, homo))

            windowCounter += 1
            windowMin      = argv.window * windowCounter
            windowMax      = windowMin + argv.window

    return

if __name__ == '__main__':
    main()
