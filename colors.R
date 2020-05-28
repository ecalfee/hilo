library(inlmisc)
library(RColorBrewer)
library(dichromat)
# color palettes for plots (maize yellow/orange, mexicana purple, parviglumis green)
# note: maize and mexicana need color pairs light/dark
#display.brewer.all(n=NULL, type="div",
#                   colorblindFriendly=T)
#display.brewer.pal(5, "BrBG")
#display.brewer.pal(3, "RdBu")
#display.brewer.pal(3, "BrBG")
#https://www.r-bloggers.com/tol-color-schemes/
#plot(inlmisc::GetColors(11, scheme = "sunset"))
#plot(inlmisc::GetColors(9, scheme = "muted", grey = T))
col_maize_mex_parv = inlmisc::GetColors(9, scheme = "muted", grey = T)[c(3,9,4)] # or c(3,9,7)
names(col_maize_mex_parv) <- c("maize", "mexicana", "parviglumis")
#plot(inlmisc::GetColors(9, scheme = "muted", grey = T, blind = "m"))
#plot(inlmisc::GetColors(2, scheme = "sunset"))
#plot(inlmisc::GetColors(11, scheme = "sunset", blind = "m"))
col_pos_neg = inlmisc::GetColors(11, scheme = "sunset")[c(3,10)] # or 10 and 3
#col_pos_neg = brewer.pal(11, "RdBu")[c(9, 2)]
#col_pos_neg = brewer.pal(11, "RdYlBu")[c(10, 3)]
#col_pos_neg = brewer.pal(3, "BrBG")[c(3, 1)]
names(col_pos_neg) <- c("+", "-")

#plot(inlmisc::GetColors(11, scheme = "sunset")[c(2,10)])
#plot(inlmisc::GetColors(9, scheme = "BuRd"))
#plot(inlmisc::GetColors(256, scheme = "BuRd"))
#plot(inlmisc::GetColors(  9, scheme = "PRGn"))
#plot(inlmisc::GetColors(256, scheme = "PRGn"))

