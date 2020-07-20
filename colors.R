library(RColorBrewer)
#library(inlmisc)
#library(dichromat)

# color palettes for plots (maize yellow/orange, mexicana purple, parviglumis green)
# note: maize and mexicana need color pairs light/dark

#https://www.r-bloggers.com/tol-color-schemes/
# explore some colors:
#plot(inlmisc::GetColors(9, scheme = "muted", grey = T, blind = "m"))
#plot(inlmisc::GetColors(2, scheme = "sunset"))
#plot(inlmisc::GetColors(11, scheme = "sunset", blind = "m"))

#col_maize_mex_parv = inlmisc::GetColors(9, scheme = "muted", grey = T)[c(3,9,4)] # or c(3,9,7)
col_maize_mex_parv = c("#DDCC77", "#AA4499", "#117733") # colors from above
names(col_maize_mex_parv) <- c("maize", "mexicana", "parviglumis")

#col_pos_neg = inlmisc::GetColors(11, scheme = "sunset")[c(3,10)] # or 10 and 3
col_pos_neg = c("#6DA5CC", "#DC3D2D") # from above
names(col_pos_neg) <- c("+", "-")

# check colors for color-blind friendly
#dichromat(colours = col_pos_neg, "deutan") #type = c("deutan", "protan", "tritan")
