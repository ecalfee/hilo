# useful functions for recombination map
source("calcMapRate.R")
source("calcMapPosCM.R")
source("calcMapPosBP.R")
source("calcMapLinearApprox.R") # cM2bp() and bp2cM() linear interpolation functions, also loads extended recombination map

# test different rmap functions for a few positions
#calcMapRate(chr = 1 , pos = 412, rmap = rmapEXT)
#calcMapRate(chr = 1 , pos = 1000, rmap = rmapEXT)
#calcMapPosCM(chr = 2, pos_bp = 500000, rmap = rmapEXT)
#calcMapPosCM(chr = 2, pos_bp = 5000000, rmap = rmapEXT)
#calcMapPosBP(chr = 2, pos_cM = 2, round = T, rmap = rmapEXT)
#calcMapPosBP(chr = 2, pos_cM = 109, round = F, rmap = rmapEXT)
#calcMapPosBP(chr = 2, pos_cM = 8888888, round = F, rmap = rmapEXT)
#calcMapPosBP(chr = 2, pos_cM = -70, round = F, rmap = rmapEXT)