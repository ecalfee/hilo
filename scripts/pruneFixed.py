# this script prunes a set of SNPs, read from a gzipped tab deliminated file,
# so that no two are closer than L bp (to reduce LD)

# to run: hilo/scripts/$ python pruneFixed.py True 0 1 100 ../data/TEST/10SNPS.mafs.gz ../data/TEST/10SNPS.pruned.mafs.gz
import sys
import csv
import gzip

# get input 
hasHeader = eval(sys.argv[1]) # is there a header to the original SNP position file? True of False
colChr = int(sys.argv[2]) # column number of SNP chromosome/scaffold/LG. (index starts at zero!)
colN = int(sys.argv[3]) # column number of SNP position column (index starts at zero!)
minL = int(sys.argv[4]) # minimum length (in bp) between 2 SNPs that are kept 
fileIn = sys.argv[5] # file/path to read in (gzipped)
fileOut = sys.argv[6] # file/path to write results to (gzip compressed)

print("Pruning: has header=" + sys.argv[1] + " colChr=" + sys.argv[2] + " colN=" + sys.argv[3] + " minL=" + sys.argv[4] + " fileIn=" + fileIn + " fileOut=" + fileOut)

pos1 = None # starting position
chr1 = None # starting chromosome
with gzip.open(fileOut, mode = "wt") as file1:
    writer = csv.writer(file1, delimiter = "\t")
    with gzip.open(fileIn, mode = "rt") as file2:
        reader = csv.reader(file2, delimiter = "\t")
        if hasHeader:
            writer.writerow(next(reader)) # prints header to output file
        for row in reader:
            if row[colChr] == chr1 and pos1 != None and (int(row[colN]) - pos1) < minL: # if same chromosome and too close, skip SNP
                #print("skip position", row[colN])
                pass
            else:
                #print("write position", row[colN])
                chr1 = row[colChr] # update current chromosome
                pos1 = int(row[colN]) # update current position
                writer.writerow(row) # print line (SNP included) 

print("Done pruning.")


