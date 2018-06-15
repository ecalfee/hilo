import pandas
import calcMapPos
rmapALL = pandas.read_csv("../data/linkage_map/ogut_fifthcM_map_agpv4.txt", sep = "\t", header = None, names = ["name", "marker", "pos_cM", "chrom", "pos_bp"])
# at last SNP
print(calcMapPos.calcMapPos(chrom = 10, pos = 150873537, rmap = rmapALL) == 106.4)
# just over
print(calcMapPos.calcMapPos(chrom = 10, pos = 150873538, rmap = rmapALL) > 106.4)
# just under
print(calcMapPos.calcMapPos(chrom = 10, pos = 150873536, rmap = rmapALL) < 106.4)
# within map
# in between:
# 6242    S_1862841  M6842    -2.2     10    1428910
print(calcMapPos.calcMapPos(chrom = 10, pos = 1428910, rmap = rmapALL) == -2.2)
print(calcMapPos.calcMapPos(chrom = 10, pos = (1477852+1428910)/2, rmap = rmapALL) == -2.1)
# way under:
print(calcMapPos.calcMapPos(chrom = 5, pos = 14, rmap = rmapALL) == (154.4 - 35.6) / (223716975-75378623) * (14 - 75378623) + 35.6)





