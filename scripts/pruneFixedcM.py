# this script prunes a set of SNPs, read from a gzipped tab deliminated file,
# so that no two are closer than L cM (RECOMBINATION DISTANCE -- to reduce LD)

# to run: hilo/scripts/$ python pruneFixedcM.py 0.5 55 0.1 ../data/TEST/10SNPS.prunedR.mafs ../data/TEST/10SNPS.mafs.gz 
import sys
import csv
import gzip
import pandas
import calcMapPos # helper function to calculate map position in cM from bp pos


# get input 
minL = int(sys.argv[1]) # minimum length (in cM) between 2 SNPs that are kept 
minInd = int(sys.argv[2]) # minimum number of individuals wiht data
minMAF = int(sys.arg[3]) # column number specifying 
fileOut = sys.argv[4] # file/path to write results to (.gz will be added for gzip compressed)
listIn = sys.argv[5:] # list of files/paths to read in (gzipped) 

print("Pruning R: minL= " + str(minL) + "cM, minInd= " + str(minInd) + " minMAF= " + str(minMAF) + " fileOut= " + fileOut + " listIn= " + listIn)

# set values
# chromo	position	major	minor	phat	nInd
colChr = 0 # column number of SNP chromosome/scaffold/LG. (index starts at zero!)
colPosBP = 1 # column number of SNP position column (index starts at zero!)
colMAF = 4 # column number specifying minimum allele freq. of SNPs to include
colMinInd = 5 # column number specifying number of individuals with data to include a SNP

# get linkage map
rmapALL = pandas.read_csv("../data/linkage_map/ogut_fifthcM_map_agpv4.txt", sep = "\t", header =    None, names = ["name", "marker", "pos_cM", "chrom", "pos_bp"])

pos1 = None # starting position
chr1 = None # starting chromosome
# first output file has all included SNPs with original bp positions, and is gzipped
with gzip.open(fileOut + ".gz", mode = "wt") as file1:
    writerSNP = writer(file1)
    # second output file has recombination position only for all included SNPs
    with open(fileOut + ".rpos", mode = "w") as file2:
        writerPOS = csv.writer(file2, delimiter = "\t")
        # for each input file, find SNPs to include
        for fileIn in listIn:
            with gzip.open(fileIn, mode = "rt") as fileRead:
                reader = csv.reader(fileRead, delimiter = "\t")
                for row in reader:
                    try (pos_bp2 = int(row[colPos])):
                        if row[colMAF] < minMAF or row[colMinInd] < minInd:
                            pass # doesn't meet minimum criteria
                        else:
                            # get new chrom
                            chr2 = row[colChr] 
                            # get new recombination position 
                            pos2 = calcMapPos.calcMapPos(chrom = chr2, pos = row[colPosBP], rmap = rmapALL)
                            if chr2 == chr1 and (pos2 - pos1) < minL:
                                if pos1 != None & pos2 < pos1:
                                    raise ValueError("Positions must be in order on each chromosome! Error at Chr" + chr2 + ":" + pos2)
                                pass # if same chromosome and too close, skip SNP
                            else: # include if SNP passes minimum criteria and not too close
                                writerSNP.writerow(row) # print line (SNP included)
                                if chr2 != chr1: # first SNP included from a chromosome
                                    writerPOS.write(1) # arbitrary but starts with difference = 1
                                else: # print difference in Morgans
                                    writerPOS.write((pos2 - pos1)/100) 
                                chr1 = chr2 # update current chromosome
                                pos1 = pos2 # update current position
 
                    except ValueError: # if row[colPos] isn't an integer, skip and go to next line because current line is a header
                        print("skipping header line: " + row[colPos])
                        continue
            
print("Done pruning.")


