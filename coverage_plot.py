'''
This is made so that we can look at genome wide coverage more easily for R plotting.
Should output 2 files:

file 1 fields:
+ scaffold
+ bin midpoint
+ overflow midpoint (for plotting scaffolds next to each other)
+ mean coverage of bin

file 2 should have the follow fields:
+ scaffold label (of everything after)
+ breakpoint = breakpoints between chromosomes

The bedfiles must account for EACH site, otherwise the coverage data will be wonky
i.e. use the -bga option
'''
###############
### Modules ###
###############
import argparse


#################
### Arguments ###
#################
parser = argparse.ArgumentParser()

parser.add_argument('-w', '--window', default = 10000, type = int, help = 'window size')
parser.add_argument('-b', '--bed', help	= 'bed file with coverage data')
parser.add_argument('-o', '--out', help	= 'outfile prefix')

argv = parser.parse_args()


#################
### Functions ###
#################
#def scaf_max_size(fai):
#	'''
#	Input
#		fai file or 2 column file with scaffold and scaffold size
#
#	Output
#		dictionary with key,value = scaffold, size
#	'''
#	scafSizeDict = {}
#	
#	with open(fai, 'r') as faiFile:
#		for i in faiFile:
#			line = i.strip().split()
#			scafSizeDict[line[0]] = int(line[1])
#
#	return scafSizeDict

def main():
    # read in data files
    bedFile = open(argv.bed, 'r')
    covFile = open('{}_{}_cov.txt'.format(argv.out, argv.window), 'w')
    brkFile = open('{}_{}_brk.txt'.format(argv.out, argv.window), 'w')

    # generate scaffold size dictionary
    # scafSizeDict = scaf_max_size(argv.fai) < might not need, since I base on scaffold

    # generate coverage dict to fill up
    # coverage dictionary = {bin #:[scaffold, coverage, bin size, mid position, overflow mid position]}
    coverageDict = {}

    # initialize overflow variable to add to lengths of subsequent contigs and also scaffold
    overflow = 0 # for when we have a new scaffold
    prevScaf = ''
    prevEnd  = 0
    #wdwflow  = 0
    # wdw flow is there to make a new bin when we're looking at a new scaffold, since we're using flow control purely
    # on position and the floor function, regardless of scaffold

    for b in bedFile:
        line  = b.strip().split()
        scaf  = line[0]
        start = int(line[1])
        end   = int(line[2]) # non-inclusive       
        cov   = int(float(line[3]))

        # give first value to scaffold and set first breakpoint
        if prevScaf == '':
            prevScaf = scaf
            brkFile.write('{}\t{}\n'.format(scaf, overflow))

        # modify overflow as needed
        if prevScaf != scaf:
            overflow += prevEnd
            prevScaf = scaf
            #wdwflow  += 1
            brkFile.write('{}\t{}\n'.format(scaf, overflow))
        else:
            prevEnd = end

        for c in range(start, end):
            window = scaf + str(c//argv.window)
            if window not in coverageDict:
                coverageDict[window] = [scaf,0,0,0,0]

            currentMidPoint = (c + 1 + (c//argv.window * argv.window))/2

            coverageDict[window][1] += cov
            coverageDict[window][2] += 1
            coverageDict[window][3] = currentMidPoint
            coverageDict[window][4] = currentMidPoint + overflow

    for j in coverageDict:
        scaf             = coverageDict[j][0]
        meanCoverage     = coverageDict[j][1]/coverageDict[j][2]
        midPointGeneral  = coverageDict[j][3]
        midPointOverflow = coverageDict[j][4]
        covFile.write('{}\t{}\t{}\t{}\n'.format(scaf, meanCoverage, int(midPointGeneral), int(midPointOverflow)))

    
    return 'BLEGH'

if __name__ == '__main__':
	main()
