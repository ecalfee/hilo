# this script calculates the map position, in cM, for the start end and midpoint
# of each CoDing Sequence (CDS) in the maize genome

# to run: hilo/scripts$ python3 getMapPosCDS.py
import csv
import pandas
import calcMapPos # helper function to calculate map position in cM from bp pos


# get input file of CDS and their bp positions
nameIn = "../data/refMaize/geneAnnotations/CDS.txt"
nameOut = "../data/refMaize/geneAnnotations/CDS_mapPos.txt" # name output file

# get linkage map
rmapALL = pandas.read_csv("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt", sep = "\t", header = None, names = ["name", "marker", "pos_cM", "chrom", "pos_bp"])

print("starting to calculate positions of CDS")

with open(nameOut, mode = "wt") as fileOut:
    writer = csv.writer(fileOut, delimiter = "\t")
    writer.writerow(["chrom", "mid_bp", "mid_cM", "length_bp", "CDS", "start_bp", "start_cM", "end_bp", "end_cM"])
    # for each input file, find SNPs to include
    with open(fileIn, mode = "rt") as fileRead:
        reader = csv.reader(fileRead, delimiter = "\t")
        for row in reader:
            chrom = row[0]
            CDS = row[1]
            start_bp = int(row[2])
            end_bp = int(row[3])
            # get midpoint of CDS and length
            mid_bp = (float(start_bp) + float(end_bp))/2.0 # may not be whole bp
            length_bp = int(end_bp) - int(start_bp)

            # get recombination positions
            start_cM = calcMapPos.calcMapPos(chrom = chrom, pos = start_bp, rmap = rmapALL)
            end_cM = calcMapPos.calcMapPos(chrom = chrom, pos = end_bp, rmap = rmapALL)
            mid_cM = calcMapPos.calcMapPos(chrom = chrom, pos = mid_bp, rmap = rmapALL)
            #print("at CDS " + str(chrom)+ ":" + str(mid_bp) + " at "+ str(mid_cM) + "cM")

            # print to file
            writer.writerow([chrom, str(mid_bp), str(mid_cM), str(length_bp), CDS, str(start_bp), str(start_cM), str(end_bp), str(end_cM)])

print("Done")
