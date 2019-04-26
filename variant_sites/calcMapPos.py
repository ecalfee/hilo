# this script has a helper function to calculate the absolute map position
# for any SNP, given a chromosome and bp position and a recombination map rmap # with bp and cM positions for some fixed width positions (e.g. Ogut .2cM map)
import pandas
# to use in other scripts: import calcMapPos

def calcMapPos(chrom, pos, rmap):
    # subset recombination map to only relevant chromosome
    rmapChr = rmap[rmap["chrom"] == chrom]
    # if position is outside the range of the map 
    # (less than 1st value or greater than last value)
    # use avg. recombination rate for the closest bin (at that end of the chr)
    if (pos < rmapChr["pos_bp"].iloc[0]):
        #print("off map left")
        left = rmapChr.iloc[0, :] # leftmost mapped position
        right = rmapChr.iloc[1, :] # next position from furthest left
    elif (pos >= rmapChr["pos_bp"].iloc[-1]):
        left = rmapChr.iloc[-2, :] # second to last position from the right
        right = rmapChr.iloc[-1, :] # rightmost mapped position on chromosome
    # otherwise if SNP is within the map, find the closest mapped positions above and below pos
    else:
        left = rmapChr[rmapChr["pos_bp"] <= pos].iloc[-1] # last mapped position smaller than or equal to pos
        right = rmapChr[rmapChr["pos_bp"] > pos].iloc[0] # first mapped position greater than pos

    # calculate cM/bp in local region between 2 closest mapped positions
    rate = (right["pos_cM"] - left["pos_cM"]) / (right["pos_bp"] - left["pos_bp"])
    #print("cM/bp rate is " + str(rate)) 
    # use this rate to find distance from leftmost position
    # note that the bp difference * rate will be negative if the position is off the map to the left of all mapped positions
    map_pos = left["pos_cM"] + (pos - left["pos_bp"])*rate 
    
    if rate <= 0: # should always be true with cleaned recombination rate map
        raise ValueError("Recomb. rate calculated is NOT strictly positive, pos = " + str(pos))

    return(map_pos) # returns position in cM (note: may be negative! b/c mapped positions go negative from some previous 0)
        


