# in this script I make a map of the sampling sites in HILO dataset for admixed populations
# packages for working with raster data
library(rgdal)
library(raster)
library(scales)
library(viridisLite) # colors
library(geosphere)
library(ade4) # mantel test for matrix correlations


# populations sampled
d = read.csv("data/PlateNames.csv", header = T, stringsAsFactors = F)
table(d$Population)
d[d$Population=="Pururandiro",] # maybe an allopatric pop (?); likely misspelling Puruandiro
tapply(d$Pop, d$Population, unique)
# gps locations for admixed/sympatric populations ** note some locations are tentative (google maps) until I get official gps coordinates
gps = read.csv("data/admixed_gps.csv", header = T, stringsAsFactors = F)

m <- read.table("data/riplasm/gps_and_elevation_for_sample_sites.txt",
                sep = "\t", header = T) %>%
  mutate(group = paste(symp_allo, zea, sep = "_")) %>%
  mutate(pop = paste0("pop", popN))

# change in plot settings needed to plot raster images
graphics.off() 
par("mar") 
par(mar=c(1,1,1,1))

# load elevational data
#mex <- raster("data/elevation_global/e10g")
mex <- raster("data/elevation_global/15-H.tif")
#my_colors = colorRampPalette(colors = c("blue", "green", "yellow", "red", "purple", "black"))(500)

#my_colors[1:6] <- "#FFFFFF"
my_colors[1:3] <- "#000666" # turn ocean dark (points close to zero elev.)
#plot(mex, col = topo.colors)
#plot(mex, col = terrain.colors)
plot(mex, col = my_colors)
#e <- drawExtent()
#class       : Extent 
#xmin        : -118.077 
#xmax        : -85.42932 
#ymin        : 13.84779 
#ymax        : 32.47407 
e <- c(-118.077, -85.42932, 13.84779, 32.47407)
SpatialPoints(coords = cbind(c(-118.077, -85.42932), c(13.8477, 32.47407)),
              proj4string = CRS("+proj=longlat +datum=WGS84"))
zoomIn <- crop(mex, e)
# plot just mexico
tiff("plots/elevation_hilo_sympatric_pts.tiff", 
       height = 2.5, width = 5, units = 'in', pointsize = 10, 
     res=1000, compression = "lzw")
par(mar=c(0,0,0,2)) # change margins
plot(zoomIn, col = my_colors)

# add GPS sampling location points
plot(SpatialPoints(coords = gps[,c("long", "lat")], 
                   proj4string = CRS("+proj=longlat +datum=WGS84")), 
     add = T, col = alpha("black", 1), pch = 20, cex = 1)
dev.off() # save


tiff("plots/elevation_hilo_sympatric_pts_viridis.tiff", 
     height = 2.5, width = 5, units = 'in', pointsize = 10, 
     res=1000, compression = "lzw")
par(mar=c(0,0,0,2)) # change margins
my_colors_viridis = viridis(n = 100, alpha = 1, begin = .1, end = .8, direction = -1, option = "D")
my_colors_viridis[1] <- "white"
plot(zoomIn, col = my_colors_viridis)
#my_colors_viridis[1] <- viridis(n = 1, alpha = 1, begin = 0, end = .2, direction = -1) # turn ocean dark

# add GPS sampling location points
#plot(SpatialPoints(coords = m[m$zea=="mexicana", c("LON", "LAT")], 
#                   proj4string = CRS("+proj=longlat +datum=WGS84")), 
#     add = T, col = ifelse(m$symp_allo[m$zea == "mexicana"] == "allopatric", "white", "black"), pch = 1, cex = 1)
plot(SpatialPoints(coords = m[m$zea=="maize", c("LON", "LAT")], 
                                      proj4string = CRS("+proj=longlat +datum=WGS84")), 
                        add = T, col = "white", pch = 6, cex = 1)
dev.off() # save. note: some points are so close together it's hard to see


# use a projection to get distance 'as the crow flies'
# between two points
distm(m[m$zea=="maize", c("LON", "LAT")][1, ], 
      m[m$zea=="maize", c("LON", "LAT")][3, ], 
      fun = distVincentyEllipsoid) # very accurate; slighly slower than spherical methods, e.g. distHaversine
# get matrices for the geographic and elevational differences
# between populations
# goal = compare geographic and elevational distance matrices to K matrix
maize_pops_byElev <- m[m$zea == "maize", ] %>% # note: may not be in same order as K matrix -- will need to reorder if so
  arrange(ELEVATION)
distmatrix <- matrix(0, 14, 14) # in meters
distelev <- matrix(0, 14, 14)
for (i in 1:nrow(maize_pops)){
  for (j in 1:nrow(maize_pops)){
    distmatrix[i,j] <- distm(maize_pops[i, c("LON", "LAT")], 
                             maize_pops[j, c("LON", "LAT")], 
                             fun = distVincentyEllipsoid)
    distelev[i,j] <- abs(maize_pops[i, "ELEVATION"] - 
                           maize_pops[j, "ELEVATION"])
  }
}
# quick way to get lower matrix of distances
distelev - as.matrix(dist(maize_pops$ELEVATION)) # all zeros
# note dist() returns a distance object, lower triangle, but as.matrix converts to full matrix

# slight positive correlation
mean((distmatrix-mean(distmatrix))*(distelev - mean(distelev)))/
  sqrt(mean((distmatrix-mean(distmatrix))^2)*
              mean((distelev - mean(distelev))^2))
# mantel permutation test to see if the two matrices are correlated
mantel.test(distmatrix, distelev, nrepeat = 100)

# alternative colors
my_colors = colorRampPalette(colors = c("blue", "green", "yellow"))(500)
my_colors[1:3] <- "#000666" # turn ocean dark
plot(zoomIn, col = my_colors)
plot(SpatialPoints(coords = gps[,c("long", "lat")], 
                   proj4string = CRS("+proj=longlat +datum=WGS84")), 
     add = T, col = alpha("black", 1), pch = 20, cex = 1)
