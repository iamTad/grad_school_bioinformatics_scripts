# Bioinformatics Scripts
These scripts are here to help me with some day to day stuff. Generally, you can figure out what they do by looking at the commments I make in the head of each script. You can look at some options with **python <script> -h** Here's a quick rundown of what each script does:

## General use
**coverage_plot.py** takes in a bedgraph file of coverage (or really anything that assigns a quantity to a linear region) and spits out mean coverage for non-overlapping sliding windows. In case you want to concatenate contigs/scaffolds/chromosomes for ease of plotting in R, the last column automatically does so (compare the last and second to last columns when you switch scaffolds).

**repeat_along_genome.py** will give you the average repeat content along the genome in non-overlapping sliding windows. Can technically calculate single nucleotide content if you put that down for "missing". Additionally, this script is case sensitive, so make sure the character you're searching for are all upper or lower case.

**reverse_complement.py** will give you a fasta file with sequences that you want to get the reverse complement of. You have to give it an initial fasta file and a list of sequences to modify.

**filter_good_mRNA.py** will filter out mRNA sequences that either doesn't have a start codon, has a premature stop codon, or doesn't have a stop codon. Does not check reverse complement, because I wrote this specifically for MAKER outputs and MAKER already tries to make all sequences start with an ATG.

## Population genetics
**hka_test.py** takes in a vcf file and calculates hka between two species in non-overlapping sliding windows across the genome. Will also perform a homogeneity test for the two species of interest. I wrote this because I couldn't find an HKA test calculator online that processes a vcf file. I used PANDAS for convenience, but that means the memory usage will be enormous--I might rewrite this using only hash tables.

The way it works, generally is:
W = number of polymorphisms in species A
X = number of polymorphisms in species B
Y = number of fixed differences between species A and species B+C (outgroup)
Z = number of fixed differences between species B and species A+C

calculate each of the above per window and also get total values for each of the above.

HKA(speciesA) = -log10(pval(chi2([W,Y],[WGlobal,ZGlobal])))
homogeneity(speciesA, speciesB) = -log10(pval(chi2([W,Y],[X,Z])))
