# this script prunes a set of SNPs, read from a tab deliminated variant sites file,
# so that no two are closer than L cM (RECOMBINATION DISTANCE -- to reduce LD)

# to run: hilo/scripts/$ python3 pruneFixedcM.py 0.0001 out_file in_file1 in_file2
import sys
import csv
import pandas
import calcMapPos # helper function to calculate map position in cM from bp pos


# get input
minL = float(sys.argv[1]) # minimum length (in cM) between 2 SNPs that are kept
fileOut = sys.argv[2] # file/path to write results to as a new sites file .var.sites and .distM distance in Morgans file
listIn = sys.argv[3:] # list of files/paths to read in (sites files, not gzipped)

print("Pruning R: minL= " + str(minL) + "cM, " + " fileOut= " + fileOut + " listIn= " + str(listIn))

# sites input file should have the following format:
# chromo	position	major	minor
colChr = 0 # column number of SNP chromosome/scaffold/LG. (index starts at zero!)
colPosBP = 1 # column number of SNP position column (index starts at zero!)

# get linkage map
rmapALL = pandas.read_csv("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt", sep = "\t", header = None, names = ["name", "marker", "pos_cM", "chrom", "pos_bp"])

pos1 = None # starting position
chr1 = None # starting chromosome
# first output file has all included SNPs with original bp positions, and is gzipped
with open(fileOut + ".var.sites", mode = "wt") as fileSNP:
    writerSNP = csv.writer(fileSNP, delimiter = "\t")
    # second output file has distance in Morgans between all included SNPs
    # (where first SNP of any chromosome is set to 1, an arbitrary choice from ancestry_hmm)
    with open(fileOut + ".distM", mode = "w") as writerPOS:
        # for each input file, find SNPs to include
        for fileIn in listIn:
            with open(fileIn, mode = "rt") as fileRead:
                reader = csv.reader(fileRead, delimiter = "\t")
                for row in reader:
                    try:
                        posBP2 = int(row[colPosBP])
                    except ValueError: # if row[colPosBP] isn't an integer, skip and go to next line because current line is a header
                        print("skipping header line: " + str(row))
                        continue

                    # get new chrom
                    chr2 = int(row[colChr])
                    # get new recombination position
                    pos2 = calcMapPos.calcMapPos(chrom = chr2, pos = posBP2, rmap = rmapALL)
                    #print("at SNP " + str(chr2)+ ":" + str(pos2) + "cM")
                    if chr2 == chr1 and (pos2 - pos1) < minL: # if same chromosome and too close, skip SNP
                        if pos1 is not None and pos2 < pos1:
                            raise ValueError("Positions must be in order on each chromosome! Error at Chr" + str(chr2) + ":" + str(pos2) + "at bp-position " + str(postBP2) + " should not come after Chr" + str(chr1) + ":" + str(pos1) + " at bp-position " + str(posBP1) + " file " + fileIn)
                        pass
                    else: # include if on diff chrom or not too close
                        writerSNP.writerow(row) # print line (SNP included)
                        if chr2 != chr1: # first SNP included from a chromosome
                            writerPOS.write(str(1) + "\n") # arbitrary but starts with difference = 1
                        else: # print difference in Morgans
                            writerPOS.write(str((pos2 - pos1)/100) + "\n")
                        chr1 = chr2 # update current chromosome
                        pos1 = pos2 # update current position (cM)
                        posBP1 = posBP2 # update current bp position
print("Done pruning by fixed cM distance.")