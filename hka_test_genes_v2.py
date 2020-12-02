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

Note:
- The vcf input is unique in that it's actually merged with a gff file using bedtools. The last column will have the gene name
- Additionally, I must manually add the vcf header to the file, something like:
> cat <(grep '#' vcf_file) modified_gff_file > input_file

- the gene ranges are...

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
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency


### Arguments ###
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', help = 'input vcf file')
parser.add_argument('-p', '--pop', help = 'list of populations')
parser.add_argument('-c', '--cut', help = 'how many missing individuals max?', type = float, default = 0.4)
parser.add_argument('-o', '--out', help = 'output file name')
parser.add_argument('-b', '--bed', help = 'is this vcf merged with gff or bed? args are bed/gff', type = str, default = 'gff')
argv = parser.parse_args()


### Functions ###
def parse_populations(popFile):
    '''
    input is the population file
    output a dictionary of individual:species
    '''
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


def get_counts(vcf, pop_dict, pop_list):
    '''
    Inputs:
    vcf = vcf file
    pop_dict = dictionary of populations:sample
    pop_list = list of populations. First two elements are sister species, third is outgroup

    Outputs:
    a dictionary of gene: {w:, x:, y:, z:}
    should also include total[w,x,y,z]
    '''
    ### initialize output dictionary
    gene_dict = {}
    # add the totals entry into the dictionary
    gene_dict[('total', 0, 0, 0)] = {'w':0, 'x':0, 'y':0, 'z':0}

    with open(vcf, 'r') as infile:
        for i in infile:
            if i.startswith('##'):
                continue
            elif i.startswith('#'):
                header = i.strip().split()
                continue
            else:
                pass

            ### read the line
            line = i.strip().split()

            ### get gene name and start/end
            gene  = line[-1].split('=')[-1] # the gene name is the last column and has an = sign before it
            scaf  = line[0]
            if argv.bed == 'gff':
                start = int(float(line[len(header) + 3]))
                end   = int(float(line[len(header) + 4]))
            elif argv.bed == 'bed':
                start = int(float(line[len(header) + 1]))
                end   = int(float(line[len(header) + 2]))

            ### add the gene to the dictionary
            if (gene, scaf, start, end) not in gene_dict:
                gene_dict[(gene, scaf, start, end)] = {'w':0, 'x':0, 'y':0, 'z':0}

            ### now start adding up counts for that gene
            ## get values for each population
            # genotypes will alwways be F/G, I can just parse by indexing the first and last element
            pop_A_allele = ''.join([line[header.index(a)][0] + line[header.index(a)][2] for a in pop_dict[pop_list[0]]]) # eg Am
            pop_B_allele = ''.join([line[header.index(a)][0] + line[header.index(a)][2] for a in pop_dict[pop_list[1]]]) # eg Nt
            pop_C_allele = ''.join([line[header.index(a)][0] + line[header.index(a)][2] for a in pop_dict[pop_list[2]]]) # eg Kl

            print(gene, pop_A_allele, pop_A_allele.count('.')/len(pop_A_allele), pop_B_allele, pop_B_allele.count('.')/len(pop_B_allele), pop_C_allele, pop_C_allele.count('.')/len(pop_C_allele))            

            # skip sites if the number of missing (".") values are equalto or beyond the threshold
            if (pop_A_allele.count('.')/len(pop_A_allele) >= argv.cut) \
                or (pop_B_allele.count('.')/len(pop_B_allele) >= argv.cut) \
                or (pop_C_allele.count('.')/len(pop_C_allele) >= argv.cut):
                # we dont' have enough data, so skip
                print(gene, pop_A_allele, pop_A_allele.count('.')/len(pop_A_allele), pop_B_allele, pop_B_allele.count('.')/len(pop_B_allele), pop_C_allele, pop_C_allele.count('.')/len(pop_C_allele))
                continue

            ### convert list of genotypes to set
            pop_A_allele_set = set(pop_A_allele)
            pop_B_allele_set = set(pop_B_allele)
            pop_C_allele_set = set(pop_C_allele)
            
            ### discard missing values in sets
            pop_A_allele_set.discard('.')
            pop_B_allele_set.discard('.')
            pop_C_allele_set.discard('.')
            print(pop_A_allele_set, len(pop_A_allele_set), pop_B_allele_set, len(pop_B_allele_set), pop_C_allele_set, len(pop_C_allele_set))
            ### start addint counts to the dictionary
            ## recall w =  polymorphic in species A
            ##        x =  polymorphic in species B
            ##        y =  all fixed A != B == C
            ##        z =  all fixed B != A == C
            if len(pop_A_allele_set) > 1:
                gene_dict[(gene, scaf, start, end)]['w'] += 1
                gene_dict[('total', 0, 0, 0)]['w'] += 1
            
            if len(pop_B_allele_set) > 1:
                gene_dict[(gene, scaf, start, end)]['x'] += 1
                gene_dict[('total', 0, 0, 0)]['x'] += 1

            if ((len(pop_A_allele_set) == 1) and
                (len(pop_B_allele_set) == 1) and
                (pop_A_allele_set != pop_B_allele_set)):
                # we now know A and B are fixed and not equal
                if pop_A_allele_set == pop_C_allele_set:
                    gene_dict[(gene, scaf, start, end)]['z'] += 1
                    gene_dict[('total', 0, 0, 0)]['z'] += 1
                elif pop_B_allele_set == pop_C_allele_set:
                    gene_dict[(gene, scaf, start, end)]['y'] += 1
                    gene_dict[('total', 0, 0, 0)]['y'] += 1
                else:
                    pass

    return gene_dict


def main():
    ### initialize output file
    out_file = open(argv.out, 'w')
    
    ### write header to output
    out_file.write('#gene\tchrom\tstart\tend\tw\tx\ty\tz\thka_A\thka_B\thomo\n')
    ### get the sample dictionaries
    pop_dict, pop_list = parse_populations(argv.pop)

    ### now get the counts for w,x,y,z
    gene_dict = get_counts(argv.vcf, pop_dict, pop_list)

    ### loop through the gene_dict with counts and calculate hka and homogeneity
    t = ('total', 0, 0, 0) # because I'm too lazy to type this out later
    for g in gene_dict:
        gene  = g[0]
        scaf  = g[1]
        start = g[2]
        end   = g[3]
        if gene == 'total':
            continue
        
        #print(gene_dict[g]['w'], gene_dict[g]['y'])
        ## get HKA for pop A
        if gene_dict[g]['w'] + gene_dict[g]['y'] > 0:
            #print(gene_dict[g]['w'], gene_dict[g]['x'], gene_dict[g]['y'], gene_dict[g]['z'],
            #        gene_dict[t]['w'], gene_dict[t]['x'], gene_dict[t]['y'], gene_dict[t]['z'])
            Ahka = -np.log10(chi2_contingency(np.array([[gene_dict[g]['w'], gene_dict[g]['y']], [gene_dict[t]['w'], gene_dict[t]['y']]]))[1]) 
        
        else:
            Ahka = np.nan

        ## get HKA for pop B
        if gene_dict[g]['x'] + gene_dict[g]['z'] > 0:
            Bhka = -np.log10(chi2_contingency(np.array([[gene_dict[g]['x'], gene_dict[g]['z']], [gene_dict[t]['x'], gene_dict[t]['z']]]))[1])

        else:
            Bhka = np.nan

        ## perform homogeneity test
        if (gene_dict[g]['w'] + gene_dict[g]['y'] > 0) and \
            (gene_dict[g]['w'] + gene_dict[g]['x'] > 0) and \
            (gene_dict[g]['x'] + gene_dict[g]['z'] > 0) and \
            (gene_dict[g]['y'] + gene_dict[g]['z'] > 0):
            homo = -np.log10(chi2_contingency(np.array([[gene_dict[g]['w'], gene_dict[g]['y']], [gene_dict[g]['x'], gene_dict[g]['z']]]))[1])

        else:
            homo = np.nan

        out_file.write(f'{gene}\t{scaf}\t{start}\t{end}\t{gene_dict[g]["w"]}\t{gene_dict[g]["x"]}\t{gene_dict[g]["y"]}\t{gene_dict[g]["z"]}\t{Ahka}\t{Bhka}\t{homo}\n')

    return

if __name__ == '__main__':
    main()
