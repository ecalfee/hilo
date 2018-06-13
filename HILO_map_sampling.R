# in this script I make a map of the sampling sites in HILO dataset for admixed populations
# packages for working with raster data
library(rgdal)
library(raster)
library(scales)


# populations sampled
d = read.csv("data/PlateNames.csv", header = T, stringsAsFactors = F)
table(d$Population)
d[d$Population=="Pururandiro",] # maybe an allopatric pop (?); likely misspelling Puruandiro
tapply(d$Pop, d$Population, unique)
# gps locations for admixed/sympatric populations ** note some locations are tentative (google maps) until I get official gps coordinates
gps = read.csv("data/admixed_gps.csv", header = T, stringsAsFactors = F)

# change in plot settings needed to plot raster images
graphics.off() 
par("mar") 
par(mar=c(1,1,1,1))

# load elevational data
mex <- raster("data/elevation_global/e10g")
mex <- raster("data/elevation_global/15-H.tif")
my_colors = colorRampPalette(colors = c("blue", "green", "yellow", "red", "purple", "black"))(500)

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


# alternative colors
my_colors = colorRampPalette(colors = c("blue", "green", "yellow"))(500)
my_colors[1:3] <- "#000666" # turn ocean dark
plot(zoomIn, col = my_colors)
plot(SpatialPoints(coords = gps[,c("long", "lat")], 
                   proj4string = CRS("+proj=longlat +datum=WGS84")), 
     add = T, col = alpha("black", 1), pch = 20, cex = 1)
