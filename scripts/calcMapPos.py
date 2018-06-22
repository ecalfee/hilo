# this script has a helper function to calculate the absolute map position
# for any SNP, given a chromosome and bp position and a recombination map rmap # with bp and cM positions for some fixed width positions (e.g. Ogut .2cM map)
import pandas
# to use in other scripts: import calcMapPos

def calcMapPos(chrom, pos, rmap):
    # subset recombination map to only relevant chromosome
    rmapChr = rmap[rmap["chrom"] == chrom]
    # if position is outside the range of the map 
    # (less than 1st value or greater than last value)
    # use avg. recombination rate over whole chromosome,
    # starting at leftmost position on chromosome (= 1st value)
    if (pos < rmapChr["pos_bp"].iloc[0] or pos >= rmapChr["pos_bp"].iloc[-1]):
        #print("off map")
        left = rmapChr.iloc[0, :]
        right = rmapChr.iloc[-1, :]
    # otherwise find the closest mapped positions above and below pos
    else:
        left = rmapChr[rmapChr["pos_bp"] <= pos].iloc[-1] # last mapped position smaller than or equal to pos
        right = rmapChr[rmapChr["pos_bp"] > pos].iloc[0] # first mapped position greater than pos

    # calculate cM/bp in local region between 2 closest mapped positions
    # (or whole chromosome if outside range of map)
    rate = (right["pos_cM"] - left["pos_cM"]) / (right["pos_bp"] - left["pos_bp"])
    #print("cM/bp rate is " + str(rate)) 
    map_pos = left["pos_cM"] + (pos - left["pos_bp"])*rate
    if rate <= 0:
        raise ValueError("Recomb. rate calculated is NOT strictly positive, pos = " + str(pos))

    return(map_pos) # returns position in cM (note: may be negative! b/c mapped positions go negative from some previous 0)
        


