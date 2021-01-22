#!/usr/bin/env Rscript
# working directory is hilo/
# plot map of sampling locations
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)
library(readxl)
library(cowplot)
library(scales)


# load variables from snakemake
# get output plot filenames
png_map = snakemake@output[["png_map"]]
# png_map = "map/plots/mexico_map_elev.png"

# get colors for plot
source(snakemake@input[["colors"]])
# source("colors.R")

# get sampling location metadata
load(snakemake@input[["meta"]])
# load("samples/HILO_MAIZE55_meta.RData")

# get teosinte occurrence data
teosinte_excel = snakemake@input[["teo"]]
# teosinte_excel = "data/zea_occurence/teocintles_historico_jsgetal.xlsx"

# read in teosinte occurence data (parviglumis & mexicana)
teo <- readxl::read_xlsx(teosinte_excel, 
                         sheet = "RegistrosFinal") %>%
  filter(Taxa %in% c("Zea mays parviglumis", "Zea mays mexicana")) %>%
  rename(ELEVATION = Alt,
         zea = Subespecie,
         year = Fecha) %>%
  arrange(ELEVATION) %>%
  mutate(latlonelev = paste(Latitud, Longitud, ELEVATION),
         latlon = paste(Latitud, Longitud),
         latlonelevyear = paste(Latitud, Longitud, ELEVATION, year)) %>%
  filter(!duplicated(latlonelevyear)) %>% # only keep unique occurence observations per location (lat/lon + elevation) and year
  dplyr::select(zea, ELEVATION, Estado, Latitud, Longitud, year)

# plot map of sampling locations
mexico <- map_data("world")  %>% # world data, country outlines
  filter(., region == "Mexico")

# function to scale elevation of pops to a bar across lat
bar_x = c(-115, -90)
bar_y = c(15, 32.5)
bar_teo = c(0, 3050)
get_scaled_xy_pos <- function(elev, bar_xy, bar_teo){
  length_xy = bar_xy[2] - bar_xy[1]
  length_teo = bar_teo[2] - bar_teo[1]
  pos = bar_xy[1] + elev*length_xy/length_teo
  return(pos)
}

pop_to_map <- filter(meta, group == "sympatric_maize") %>%
  arrange(desc(ELEVATION)) %>%
  mutate(LOCALITY_ELEVATION = factor(paste0(LOCALITY, " (", ELEVATION, "m)"), 
                                     ordered = T, 
                                     levels = unique(paste0(.$LOCALITY, " (", .$ELEVATION, "m)"))),
         pos_x = get_scaled_xy_pos(elev = ELEVATION,
                                   bar_xy = bar_x,
                                   bar_teo = bar_teo),
         pos_y = get_scaled_xy_pos(elev = ELEVATION,
                                   bar_xy = bar_y,
                                   bar_teo = bar_teo))

# make map
p_map <-  pop_to_map %>%
  ggplot(data = .) +
  geom_polygon(data = mexico,
    aes(x = long, y = lat, group = group), 
               fill = "lightgrey") +
  #coord_cartesian() +
  coord_equal(#xlim = c(-110, -90)
    xlim = bar_x,
    ylim = bar_y,
                ) +
  theme_classic() +
  geom_point(
    aes(x = LON,
        y = LAT,
        shape = LOCALITY_ELEVATION,
        color = LOCALITY_ELEVATION),
             size = 2,
             alpha = 1) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(direction = 1, option = "viridis") +
  labs(color = NULL, shape = NULL) +
  #scale_shape_manual(values = rep(c(15, 17, 19, 0, 2, 1), 3)) +
  #scale_shape_manual(values = rep(c(1, 17, 0, 19, 2, 15), 3))
  #scale_shape_manual(values = rep(c(1, 0, 2, 4), 5))
  scale_shape_manual(values = rep(c(2, 0, 3, 1, 6, 5), 5))
p_map

ggsave(filename = png_map,
       plot = p_map,
       device = "png", height = 4, width = 7.5, units = "in", dpi = 300)

# add scale bar and lines
p_map_lines <- pop_to_map %>%
  ggplot(data = .) +
  geom_polygon(data = mexico,
               aes(x = long, y = lat, group = group), 
               fill = "lightgrey") +
  coord_equal(
    xlim = bar_x,
    ylim = bar_y
  ) +
  scale_y_continuous(expand = c(0,0), 
                     #limits = bar_y, 
                     #labels = scales::number_format(accuracy = 0.0)
                     ) +
  scale_x_continuous(expand = c(0,0), 
                     #limits = bar_x, 
                     #labels = scales::number_format(accuracy = 0.0)
                     ) +
  theme_classic() +
  geom_point(
    aes(x = LON,
        y = LAT,
        color = ELEVATION),
    size = .1) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis(direction = -1, option = "viridis") +
  geom_segment(aes(x = LON, y = LAT, xend = pos_x,
                   color = ELEVATION), yend = bar_y[2]) +
  labs(color = "Elevation (m)")
p_map_lines

# how big do I have to make it to see all the points in mexico?
zoom_map <- pop_to_map %>%
  ggplot(data = .) +
  #geom_polygon(aes(x = long, y = lat, group = group), 
  #             fill = "lightgrey") +
  #coord_cartesian() +
  coord_equal(#xlim = c(-110, -95),
              #ylim = c(15, 28)
    xlim = c(-102.5, -98.2),
    ylim = c(18.8, 20.3)
    ) +
  theme_classic() +
  geom_point(aes(x = LON,
                 y = LAT,
                 shape = LOCALITY_ELEVATION,
                 color = LOCALITY_ELEVATION),
             size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  #scale_shape_manual(values = rep(c(15, 17, 19, 0, 2, 1), 3)) +
  scale_shape_manual(values = rep(c(2, 0, 3, 1, 6, 5), 5)) +
  scale_color_viridis_d(direction = 1, option = "viridis") +
  labs(shape = F, color = F)
# try adding greyscale elevation underneath?
# can you make a pop-out square for samples closer together?
zoom_map
ggsave(filename = "map/plots/zoom_cluster_sampling_locations.png",
       plot = zoom_map,
       device = "png", 
       height = 4, width = 7.5, 
       units = "in", dpi = 300)

# histogram of teosinte occurence data
p_teo_hist <- teo %>%
  ggplot(.) +
  geom_histogram(aes(x = ELEVATION, fill = zea),
                 alpha = 0.5, 
                 bins = 35,
                 position = "identity") +
  theme_classic() +
  scale_x_continuous(position = "top", limits = bar_teo, expand = c(0,0)) +
  scale_y_continuous(#expand = c(0,0), 
                     #limits = c(0, 125), 
                     breaks = c(0, 30, 60, 90)) +
  scale_fill_manual(values = col_maize_mex_parv) +
  ggtitle("Teosinte occurence data 1842–2016") +
  labs(fill = "Teosinte", x = "Elevation (m)", y = "Observations")
# plot(p_teo_hist)
ggsave(filename = png_teo_hist,
       plot = p_teo_hist,
       device = "png", 
       height = 3, width = 5, 
       units = "in", dpi = 300)

# p_map_lines, p_teo_hist

p_combined <- grid.arrange(grobs = list(
  ggplotGrob(p_teo_hist + 
               geom_point(data = pop_to_map,
                          aes(x = ELEVATION,
                              color = ELEVATION),
                          y = 0,
                          shape = 3) +
               scale_color_viridis(direction = -1, option = "viridis") +
               theme(plot.title = element_blank(),
                     legend.position = "None",
                     plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
    ),
  ggplotGrob(p_map_lines + 
               theme(legend.position = "None",
                     plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"))),
  cowplot::get_legend(p_teo_hist + theme(legend.position = "left")), # + labs(fill = "Historical Occurrence")),
  cowplot::get_legend(p_map + labs(color = "Sympatric sites", #element_blank(),
                                   shape = "Sympatric sites") + #element_blank()) + # Sympatric sites?
                        guides(color = guide_legend(override.aes = list(shape = 16)))) 
  ),
  layout_matrix = rbind(
    c(1,3),
    c(2,4)),
  heights = c(2, 5),
  widths = c(5, 2))

p_combined
png_map_and_teo = "map/plots/mexico_lines_elev_teo2.png"
ggsave(png_map_and_teo, 
       plot = p_combined, 
       device = "png", 
       width = 7,#5, 
       height = 5,#5, 
       units = "in",
       dpi = 300)

p_combined2 = cowplot::plot_grid(ggplotGrob(p_teo_hist + 
                                geom_point(data = pop_to_map,
                                           aes(x = ELEVATION,
                                               color = ELEVATION),
                                           y = 0,
                                           shape = 3) +
                                scale_color_viridis(direction = -1, option = "viridis") +
                                theme(plot.title = element_blank(),
                                      legend.position = "None",
                                      plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
), 
ggplotGrob(p_map_lines + 
             theme(legend.position = "None",
                   plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"))), 
                   
ncol = 1, 
align = "v")
p_combined2
p_combined
png_map_and_teo = "map/plots/mexico_lines_elev_teo.png"
ggsave(png_map_and_teo, 
       plot = p_combined, 
       device = "png", 
       width = 5, 
       height = 5, 
       units = "in",
       dpi = 300)


# ---------- geographic distance -------------- #
# get geographic distances between populations!
all_places <- meta %>%
  filter(group == "sympatric_maize") %>%
  arrange(., desc(ELEVATION)) %>%
  .$LOCALITY
edges_geo <- data.frame(from = unlist(lapply(1:length(all_places), 
                                             function(i) 
                                               rep(all_places[i], (length(all_places)-i)))),
                        to = unlist(lapply(2:length(all_places), 
                                           function(i) 
                                             all_places[i:length(all_places)])),
                        stringsAsFactors = F) %>% # list all unique pairs of locations
  left_join(., 
            dplyr::select(meta_pops_list[["maize"]], LOCALITY, LAT, LON), 
            by = c("from"="LOCALITY")) %>%
  left_join(., 
            dplyr::select(meta_pops_list[["maize"]], LOCALITY, LAT, LON), 
            by = c("to"="LOCALITY"),
            suffix = c(".from", ".to")) %>%
  # calculate pairwise geodesic distance between locations (in m)
  dplyr::mutate(., distance = geodist(x = dplyr::select(., LON.from, LAT.from),
                                      y = dplyr::select(., LON.to, LAT.to),
                                      paired = T,
                                      measure = "geodesic")) %>%
  dplyr::mutate(distance_km = distance/1000) %>%
  dplyr::mutate(distance_km_truncated = ifelse(distance_km > 500, 500, distance_km)) %>%
  dplyr::select(from, to, distance_km, distance_km_truncated)

net_tidy_geo <- tbl_graph(nodes = nodes, 
                          edges = edges_geo,
                          directed = T)  
# complement (upside down) of the linear network plot:
p_geo_dist <- ggraph(net_tidy_geo, layout = "linear") +
  geom_edge_arc(aes(width = log(distance_km)), alpha = 0.5,
                force_flip = T) +
  geom_node_point(aes(color = ELEVATION), size = 3) + # why can't I do x = ELEVATION?
  theme_graph(base_family = 'Helvetica') +
  scale_edge_width(range = c(2, 0.01)) +
  scale_color_viridis(direction = -1) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(c(t = 20, r = 5.5, b = 100, l = 5.5), 
                             unit = "pt")) +
  guides(color = F) +
  ggtitle("Geographic distance between populations")


library(mapdata)
map("worldHires","Mexico")
x <- dplyr::select(meta, LOCALITY, ELEVATION, LAT, LON, group) %>%
  filter(group == "sympatric_maize") %>%
  unique()
points(x$LON,x$LAT,pch=19)
abline(h=15)
ey<-rep(15,14)
ex=-110+(x$ELEVATION-1500)/100
points(ex,ey,pch=3)
for(i in 1:14){
  segments(ex[i],ey[i],x$LON[i],x$LAT[i])
}
text(ex+1,ey,labels=x$ELEVATION)

plot(ex,ey,pch=3)
