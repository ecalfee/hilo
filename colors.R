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
col_group_zea = c("darkgoldenrod", "#DDCC77", "maroon", "#AA4499", "#117733") # colors from above

shape_group_zea = c(1, 19, 2, 17, 0)
alphas_group_zea = c(1, .45, 1, .45, 1)
names(col_group_zea) <- c("allopatric_maize", "sympatric_maize",
                          "allopatric_mexicana", "sympatric_mexicana",
                          "parviglumis")
names(shape_group_zea) <- c("allopatric_maize", "sympatric_maize",
                            "allopatric_mexicana", "sympatric_mexicana",
                            "parviglumis")
names(alphas_group_zea) <- c("allopatric_maize", "sympatric_maize",
                             "allopatric_mexicana", "sympatric_mexicana",
                             "parviglumis")
zea_group_labels <- c("Allopatric maize", "Sympatric maize",
                      "Allopatric mexicana", "Sympatric mexicana",
                      "Parviglumis")
names(zea_group_labels) <- c("allopatric_maize", "sympatric_maize",
                             "allopatric_mexicana", "sympatric_mexicana",
                             "parviglumis")
