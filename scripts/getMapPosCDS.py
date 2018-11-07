# this script calculates the map position, in cM, for the start end and midpoint
# of each CoDing Sequence (CDS) in the maize genome

# to run: hilo/scripts$ python3 getMapPosCDS.py 1 "../data/refMaize/geneAnnotations/CDS_1.txt" "../data/refMaize/geneAnnotations/CDS_1_mapPos.txt"
import csv
import sys
import pandas
import calcMapPos # helper function to calculate map position in cM from bp pos


# get input file of CDS and their bp positions
chromo = int(sys.argv[1]) # e.g. 1, which chromosome
nameIn = sys.argv[2] # e.g. "../data/refMaize/geneAnnotations/CDS_1.txt"
nameOut = sys.argv[3] # e.g. "../data/refMaize/geneAnnotations/CDS_1_mapPos.txt" # name output file

# get linkage map
rmapALL = pandas.read_csv("../data/linkage_map/ogut_fifthcM_map_agpv4_INCLUDE.txt", sep = "\t", header = None, names = ["name", "marker", "pos_cM", "chrom", "pos_bp"])

rmapCHR = rmapALL[rmapALL["chrom"] == chromo]

print("starting to calculate positions of CDS on chromo " + str(chromo))

with open(nameOut, mode = "wt") as fileOut:
    writer = csv.writer(fileOut, delimiter = "\t")
    writer.writerow(["chrom", "mid_bp", "mid_cM", "width_bp", "start_bp", "start_cM", "end_bp", "end_cM"])
    # for each input file, find SNPs to include
    with open(nameIn, mode = "rt") as fileRead:
        reader = csv.reader(fileRead, delimiter = "\t")
        for row in reader:
            if (row[0] == "start"):
                print("skipping line: " + str(row))
                continue # skip header line
            else:
                pass

            start_bp = int(row[0])
            end_bp = int(row[1])
            width_bp = int(row[2])
            #print("chrom is " + str(chromo) + " start_bp is " + str(start_bp))

            # get midpoint of CDS and length
            mid_bp = (float(start_bp) + float(end_bp))/2.0 # may not be whole bp
            #print("mid_bp is " + str(mid_bp) + " and length_bp is " + str(length_bp))
            # get recombination positions
            start_cM = calcMapPos.calcMapPos(chrom = chromo, pos = start_bp, rmap = rmapCHR)
            #print("start_cM is " + str(start_cM))
            end_cM = calcMapPos.calcMapPos(chrom = chromo, pos = end_bp, rmap = rmapCHR)
            #print("end_cM is " + str(end_cM))
            mid_cM = calcMapPos.calcMapPos(chrom = chromo, pos = mid_bp, rmap = rmapCHR)
            #print("at CDS " + str(chromo)+ ":" + str(mid_bp) + " at "+ str(mid_cM) + "cM")

            # print to file
            writer.writerow([str(chromo), str(mid_bp), str(mid_cM), str(width_bp), str(start_bp), str(start_cM), str(end_bp), str(end_cM)])

print("Done")
