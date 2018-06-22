# this script prunes a set of SNPs, read from a list of gzipped tab deliminated file(s),
# so that no two SNPs are closer than L BP (PHYSICAL DISTANCE -- to reduce LD)

# to run: hilo/scripts/$ python3 pruneFixedcM.py 10 55 0.1 ../data/TEST/10SNPS.prunedBP.mafs ../data/TEST/10SNPS.mafs.gz ../data/TEST/10SNPS.B.mafs.gz
import sys
import csv
import gzip
import pandas

# get input 
minL = float(sys.argv[1]) # minimum length (in BP) between 2 SNPs that are kept 
minInd = int(sys.argv[2]) # minimum number of individuals with data
minMAF = float(sys.argv[3]) # minimum minor allele frequency
fileOut = sys.argv[4] # file/path to write results to (.gz will be added for gzip compressed)
listIn = sys.argv[5:] # list of files/paths to read in (gzipped) 

print("Pruning by BP: minL= " + str(minL) + "bp, minInd= " + str(minInd) + " minMAF= " + str(minMAF) + " fileOut= " + fileOut + " listIn= " + str(listIn))

# set values
# chromo	position	major	minor	phat	nInd
colChr = 0 # column number of SNP chromosome/scaffold/LG. (index starts at zero!)
colPosBP = 1 # column number of SNP position column (index starts at zero!)
colMAF = 4 # column number specifying minimum allele freq. of SNPs to include
colNInd = 5 # column number specifying number of individuals with data to include a SNP

pos1 = None # starting position
chr1 = None # starting chromosome
# first output file has all included SNPs with original bp positions, and is gzipped
with gzip.open(fileOut + ".gz", mode = "wt") as fileSNP:
    writerSNP = csv.writer(fileSNP, delimiter = "\t")
    # second output file has recombination position only for all included SNPs
    with open(fileOut + ".rpos", mode = "w") as writerPOS:
        # for each input file, find SNPs to include
        for fileIn in listIn:
            with gzip.open(fileIn, mode = "rt") as fileRead:
                reader = csv.reader(fileRead, delimiter = "\t")
                for row in reader:
                    try:
                        pos2 = int(row[colPosBP])
                    except ValueError: # if row[colPosBP] isn't an integer, skip and go to next line because current line is a header
                        print("skipping header line: " + str(row))
                        continue
                    if float(row[colMAF]) < minMAF or int(row[colNInd]) < minInd:
                        #print("no pass, MAF=" + row[colMAF] + ", NInd=" + row[colNInd])
                        pass # doesn't meet minimum criteria
                    else: # SNP passes minimum criteria
                        # get new chrom
                        chr2 = int(row[colChr])
                        #print("at SNP " + str(chr2)+ ":" + str(pos2) + "cM") 
                        if chr2 == chr1 and (pos2 - pos1) < minL: # if same chromosome and too close, skip SNP
                            if pos1 is not None and pos2 < pos1:
                                raise ValueError("Positions must be in order on each chromosome! Error at Chr" + str(chr2) + ":" + str(pos2) + " should not come after Chr" + str(chr1) + ":" + str(pos1) + " file " + fileIn)
                            pass 
                        else: # include if on diff chrom or not too close
                            writerSNP.writerow(row) # print line (SNP included)
                            # update for any included SNP 
                            chr1 = chr2 # update current chromosome
                            pos1 = pos2 # update current position (cM)
print("Done pruning by fixed BP distance.")


