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
library(geodist)


# load variables from snakemake
# get output plot filenames
png_map = snakemake@output[["png_map"]]
# png_map = "map/plots/mexico_map_elev.png"
png_teo_hist = snakemake@output[["png_teo_hist"]]
# png_teo_hist = "map/plots/teosinte_hist.png"
png_map_teo_color = snakemake@output[["png_map_teo_color"]]
# png_map_teo_color = "map/plots/mexico_lines_elev_teo_color.png"
png_map_teo_black = snakemake@output[["png_map_teo_black"]]
# png_map_teo_black = "map/plots/mexico_lines_elev_teo_black.png"

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
bar_teo = c(-0, 3050)
get_scaled_xy_pos <- function(elev, bar_xy, bar_teo){
  length_xy = bar_xy[2] - bar_xy[1]
  length_teo = bar_teo[2] - bar_teo[1]
  pos = bar_xy[1] + elev*length_xy/length_teo
  return(pos)
}

pop_to_map <- filter(meta, group == "sympatric_maize") %>%
  arrange(desc(ELEVATION)) %>%
  dplyr::select("LOCALITY", "ELEVATION", "LAT", "LON") %>%
  filter(!duplicated(.)) %>%
  mutate(LOCALITY_ELEVATION = factor(paste0(LOCALITY, " (", ELEVATION, "m)"), 
                                     ordered = T, 
                                     levels = unique(paste0(.$LOCALITY, " (", .$ELEVATION, "m)"))),
         pos_x = get_scaled_xy_pos(elev = ELEVATION,
                                   bar_xy = bar_x,
                                   bar_teo = bar_teo),
         pos_y = get_scaled_xy_pos(elev = ELEVATION,
                                   bar_xy = bar_y,
                                   bar_teo = bar_teo))


# make a map with scale bar and lines
p_map <- pop_to_map %>%
  ggplot(data = .) +
  geom_polygon(data = mexico,
               aes(x = long, y = lat, group = group), 
               fill = "lightgrey") +
  coord_equal(xlim = bar_x, ylim = bar_y) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  geom_point(
    aes(x = LON,
        y = LAT,
        color = LOCALITY_ELEVATION),
    size = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_color_viridis_d(direction = 1, option = "viridis") +
  guides(color = guide_legend(title = "Maize/Mexicana Sites", 
                       override.aes = list(linetype = 0, size = 2)))
#p_map

p_map_lines_color = p_map + 
  geom_segment(aes(x = LON, y = LAT, xend = pos_x,
                   color = LOCALITY_ELEVATION), yend = bar_y[2])
p_map_lines_black = p_map + 
  geom_segment(aes(x = LON, y = LAT, xend = pos_x),
                   color = "black", yend = bar_y[2]) +
  geom_point(
    aes(x = LON,
        y = LAT,
        color = LOCALITY_ELEVATION),
    size = 2)# to add on top of lines

# histogram of teosinte occurence data
p_teo_hist <- teo %>%
  ggplot(.) +
  geom_histogram(aes(x = ELEVATION, fill = zea),
                 alpha = 0.75,
                 binwidth = 100,
                 position = "identity") +
  theme_classic() +
  scale_x_continuous(position = "top", 
                     limits = bar_teo, # 2 missing rows are ok (this is for empty histogram bins from (-50,50) with zero height). 
                     # you can see this by changing limits to c(-50,3050) and the ggplot2 warning goes away 
                     expand = c(0, 0)
                     ) +
  scale_fill_manual(values = col_maize_mex_parv) +
  ggtitle("Teosinte occurence data 1842–2016") +
  labs(fill = "Teosinte", x = "Elevation (m)", y = "Observations")
  
# p_teo_hist
ggsave(filename = png_teo_hist,
       plot = p_teo_hist,
       device = "png", 
       height = 3, width = 5, 
       units = "in", dpi = 300)

# make combined map and histogram
p_combined_color <- grid.arrange(grobs = list(
  ggplotGrob(p_teo_hist + 
               geom_point(data = pop_to_map,
                          aes(x = ELEVATION,
                              color = LOCALITY_ELEVATION),
                          y = 0,
                          shape = 3) +
               geom_hline(yintercept = 0) +
               scale_color_viridis_d(direction = 1, option = "viridis") +
               theme(plot.title = element_blank(),
                     legend.position = "None",
                     plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt")) +
               geom_segment(x = 0, 
                            y = 118, 
                            xend = 3050, 
                            yend = 118,
                            color = col_maize_mex_parv[["maize"]],
                            arrow = arrow(length = unit(1, "mm"), 
                                          ends = "last",
                                          type = "closed")) +
               geom_text(x = 2700, 
                         y = 102, 
                         color = col_maize_mex_parv[["maize"]],
                         label = "maize range",
                         size = 3)
    ),
  ggplotGrob(p_map_lines_color + 
               theme(legend.position = "None",
                     plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"))),
  cowplot::get_legend(p_teo_hist + 
                        theme(legend.position = "left")),
  cowplot::get_legend(p_map),
  textGrob(label = "A", 
           just = "top",
           x = unit(0.5, "lines"), 
           y = unit(0, "lines")),
  textGrob(label = "B", 
           just = "top",
           x = unit(0.5, "lines"), 
           y = unit(0, "lines"))
  ),
  layout_matrix = rbind(
    c(5, 1,3),
    c(6, 2,4)),
  heights = c(2, 5),
  widths = c(.1, 5, 2))

p_combined_color
ggsave(png_map_and_teo_color, 
       plot = p_combined_color, 
       device = "png", 
       width = 7.1,#5, 
       height = 5,#5, 
       units = "in",
       dpi = 300)

p_combined_black <- grid.arrange(grobs = list(
  ggplotGrob(p_teo_hist + 
               geom_point(data = pop_to_map,
                          aes(x = ELEVATION),
                              color = "black",
                          y = 0,
                          shape = 3) +
               geom_hline(yintercept = 0) +
               theme(plot.title = element_blank(),
                     legend.position = "None",
                     plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
  ),
  ggplotGrob(p_map_lines_black + 
               theme(legend.position = "None",
                     plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"))),
  cowplot::get_legend(p_teo_hist + 
                        theme(legend.position = "left")),
  cowplot::get_legend(p_map) 
),
layout_matrix = rbind(
  c(1,3),
  c(2,4)),
heights = c(2, 5),
widths = c(5, 2))

p_combined_black
png_map_and_teo_black = "map/plots/mexico_lines_elev_teo_black.png"
ggsave(png_map_and_teo_black, 
       plot = p_combined_black, 
       device = "png", 
       width = 7,#5, 
       height = 5,#5, 
       units = "in",
       dpi = 300)


# ---------- geographic distance -------------- #
# plot heatmap of geographic distance between populations

# get geographic distances between populations!
all_places <- pop_to_map %>%
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
            dplyr::select(pop_to_map, LOCALITY, LAT, LON), 
            by = c("from"="LOCALITY")) %>%
  left_join(., 
            dplyr::select(pop_to_map, LOCALITY, LAT, LON), 
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

distances <- bind_rows(edges_geo,
          edges_geo %>%
            mutate(from_old = from,
                   from = to,
                   to = from_old)) %>%
  mutate(from = factor(from, ordered = T, levels = all_places[length(all_places):1]),
         to = factor(to, ordered = T, levels = all_places[length(all_places):1]))

dist_heatmap <- distances %>% # add 2nd half of pairs (reversed to <-> from)
  ggplot(aes(x = from, y = to, fill = distance_km_truncated)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis()





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
