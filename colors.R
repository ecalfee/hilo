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

# color and shapes for groups
#col_group_zea = c("#DDCC77", "#AA4499", "#DDCC77", "#AA4499") # colors from above
#col_group_zea = c("gold", "deeppink", "#DDCC77", "#AA4499") # colors from above
#col_group_zea = c("lightgoldenrod", "maroon", "#DDCC77", "#AA4499") # colors from above
col_group_zea = c("darkgoldenrod", "maroon", "#DDCC77", "#AA4499") # colors from above
#col_group_zea = c("darkgoldenrod", "deeppink", "lightgoldenrod", "#AA4499") # colors from above


shape_group_zea = c(1, 2, 19, 17)
#shape_group_zea = c(19, 17, 1, 2)
alphas_group_zea = c(1, 1, .45, .45)
names(col_group_zea) <- c("allopatric_maize", "allopatric_mexicana",
                              "sympatric_maize", "sympatric_mexicana")
names(shape_group_zea) <- c("allopatric_maize", "allopatric_mexicana",
                              "sympatric_maize", "sympatric_mexicana")
names(alphas_group_zea) <- c("allopatric_maize", "allopatric_mexicana",
                            "sympatric_maize", "sympatric_mexicana")
zea_group_labels <- c("Allopatric maize", "Allopatric mexicana",
                      "Sympatric maize", "Sympatric mexicana")



# alternative color options:
# color blind palette:
#cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# maize vs. mexicana (sympatric) vs. parviglumis
#scale_color_manual(values = cbPalette[c(2,4,8)])
# maize vs. mexicana (allopatric)
#scale_color_manual(values = cbPalette[c(7,6)])
