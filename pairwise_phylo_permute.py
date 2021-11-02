'''
The difference between this and the original is that this does a permutation test

This script will do a permutation test to determine how clustered something is.
+ first find the mean pairwise distance between TEs in a species
+ then shuffle the leaves in a tree using ete3 module
+ get mean pairwise distance again
+ get p value by dividing all instances where mean pairwise distance <= real data by # permutations


Input
+ newick file
+ file of species (one per line)

Output
+ file
    + repeat name
    + species 1 no_samples;mean_pairwise_distance;p-val;multipletest_p-val
    + species 2 same thing 
    + etc.
    + set to -1 if not enough

Parameters
+ percent cutoff for TE abundance
    + set as max(10%, 5) for default
+ number of permutations
    + set as 10000 as default
'''
import argparse
import random
from ete3 import Tree


def get_species(sample):
    '''
    sample = input sample file
    returns list of samples in the same order
    '''
    sample_list = []
    sample_dict = {}
    with open(sample, 'r') as ifile:
        for i in ifile:
            sample_list.append(i.strip())
            sample_dict[i.strip()] = []

    return sample_list, sample_dict


def phylo_leaf_name(new: str):
    '''
    This function will read in the newick string and give us a dictionary of
    species and sample name. I.e.:
        {Bm: [Bm_te_1, Bm_te_2, ...]}

    Args
        newick = newick string
    Returns
        dictionary of species specific leaf
        list of all leaves
    '''
    leaf_dict = {}
    leaf_list = []

    t = Tree(new)

    for leaf in t:
        leaf_name = str(leaf).split('-')[-1] # sample name

        ### initialize species in dictionary
        if leaf_name[:2] not in leaf_dict:
            leaf_dict[leaf_name[:2]] = []

        ### add sample to dictionary and list
        leaf_dict[leaf_name[:2]].append(leaf_name)
        leaf_list.append(leaf_name)

    return leaf_dict, leaf_list


def get_dist(new: str, leaf_l: list, leaf_d: dict, species_l: list,
             seed: int = 0, p: float = 0.1, c: int = 5):
    '''
    arguments
        new    = newick tree in string format
        leaf_l = list of all leaves
        leaf_d = dictionary of leaves per species
        seed   = seed if we're randomizing
                    + set to 0 if we want to not randomize
        p = minimum proportion of all TEs required (from argparse)
        c = minimum counts of TEs required (from argparse)
    returns
        dictionary of pairwise difference
    '''
    ### initialize output dictionary
    dist_dict = {k:0.0 for k in species_l}

    ### initialize tree structure
    phylo = Tree(new)

    ### randomize if needed
    if seed != 0:
        random.seed(seed) # set the seed
        leaf_shuf = random.sample(leaf_l, len(leaf_l)) # randomize samples

        # relabel trees based on randomized samples
        counter = 0 # keeps track of which element in the randomized sample we're looking at
        for leaf in phylo:
            leaf.name = leaf_shuf[counter] # rename
            counter += 1

    ### now calculate mean pairwise distance
    for species in species_l: # loop through each species
        if species not in leaf_d:
            dist_dict[species] = -1
            continue
        elif len(leaf_d[species]) < max(p*len(leaf_l), c):
            dist_dict[species] = -1
            continue

        tot_dist = [ (phylo&leaf_d[species][a]).get_distance(phylo&leaf_d[species][b]) 
                        for a in range(len(leaf_d[species]))
                        for b in range(a+1, len(leaf_d[species])) ]
        no_pairs = ( len(leaf_d[species]) * (len(leaf_d[species]) - 1) ) / 2 # n choose 2 = n*(n-1)/2
        mean_dist = sum(tot_dist) / no_pairs
        dist_dict[species] = mean_dist

    return dist_dict


def main(argv):
    ### initialize output file
    ofile = open(argv.out, 'w')

    ### get list of species
    species_list, species_dict = get_species(argv.spec)
    
    ### write the header
    species_header = '\t'.join(species_list)
    ofile.write(f'transposon\t{species_header}\n')

    ### get the tree
    with open(argv.new, 'r') as ifile:
        for i in ifile: # should only be one line
            newick = i.strip()

    ### get the sample list and dictionary
    te_d, te_l = phylo_leaf_name(newick)

    ### do the permutations
    for i in range(argv.perm): # loop through all the permutations
        pairwise_dist = get_dist(new = newick, leaf_l = te_l, leaf_d = te_d, species_l = species_list, seed = i)
        for j in pairwise_dist: # loop through all results/species
            species_dict[j].append(pairwise_dist[j]) # we should now 

    ### get number of tests (should be number of species analyzed)
    #print(species_dict)
    no_sample = len([i for i in species_dict if species_dict[i][0] != -1])

    ### get the p_value
    pval_dict = {}
    for s in species_dict: # loop through all results
        # recall that the first element is the value of the true data
        if species_dict[s][0] == -1: # -1 if we don't test this
            pval_dict[s] = (-1, -1)
            continue
        no_below_truth = len([i for i in species_dict[s] if i <= species_dict[s][0]])
        pval = no_below_truth / argv.perm
        pval_dict[s] = (pval, pval*no_sample) # second element is corrected pvalue fo rmultiple testing

    ### prepare output
    entry = [] # fill this up for all species
    for i in species_list:
        if i not in te_d:
            species_te_count = 0
        else:
            species_te_count = len(te_d[i])
        entry.append(f'{species_te_count};{species_dict[i][0]};{pval_dict[i][0]};{pval_dict[i][1]}')
    
    ### write results
    entry = '\t'.join(entry)
    ofile.write(f'{argv.new}\t{entry}\n')

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--new', help = 'input newick file')
    parser.add_argument('-o', '--out', help = 'output file')
    parser.add_argument('-s', '--spec', help = 'species file or list')
    parser.add_argument('-pm', '--perm', help = 'number of permutations to use, default 10k',
                        type = int, default = 10000)
    parser.add_argument('-p', '--prop', help = 'percentage of TEs from species default 0.1', 
                        type = float, default = 0.1)
    parser.add_argument('-c', '--count', help = 'minimum absolute number of TEs per species, default 5', 
                        type = int, default = 5)
    args = parser.parse_args()

    main(args)
