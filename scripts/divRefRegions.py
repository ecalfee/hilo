# this script divides a reference genome into smaller regions below the chromosome level

# to run: 
# mkdir ../data/refMaize/divide_50Mb
# hilo/scripts/$ python divRefRegions.py 50000000 ../data/refMaize/Zea_mays.AGPv4.dna.chr.fa.fai ../data/refMaize/divide_50Mb
import sys
import csv

# get input 
maxL = int(sys.argv[1]) # max length (in bp) for a region 
fileIn = sys.argv[2] # file/path to read in with total chromosome lengths
pathOut = sys.argv[3] # file/path to write results to

print("Dividing genome into regions size " + sys.argv[1] + "bp")
n = 0 # first chunk
print(pathOut + "/region_" + str(n) +".txt")
with open(fileIn, mode = "r") as fileRead:
    reader = csv.reader(fileRead, delimiter = "\t")
    for row in reader:
        chrom = row[0]
        totL = int(row[1]) # total length of chromosome
        if chrom == "Pt" or chrom == "Mt": # ignore non autosomal DNA
            continue
        else:
            start = 0
            end = start + maxL - 1
            while totL > end: # iterate until a new region spills over total length
                with open(pathOut + "/region_" + str(n) +".txt", mode = "w") as f:
                    f.write(chrom + ":" + str(start) + "-" + str(end) + "\n")
                # update region number, start position (one past previous end) 
                # and end position
                n = n + 1
                start = start + maxL
                end = start + maxL - 1
                if n > 1000:
                    print("exiting -- tried to create over 1000 files")
                    break # problem
            # when end goes over totL, write one more region for that chromosome  
            with open(pathOut + "/region_" + str(n) +".txt", mode = "w") as f:
                f.write(chrom + ":" + str(start) + "-" + str(end) + "\n")
            # and update variables
            n = n + 1
            start = None # new chromosomes start over in position
            end = None


print("Done, # total regions = " + str(n-1))


