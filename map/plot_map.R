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
  guides(color = guide_legend(title = "Maize/Mexicana sites", 
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
                     plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
    ),
  ggplotGrob(p_map_lines_color + 
               theme(legend.position = "None",
                     plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt"))),
  cowplot::get_legend(p_teo_hist + 
                        labs(fill = "Teosinte ranges") +
                        theme(legend.position = "left")),
  cowplot::get_legend(p_map),
  textGrob(label = "A", 
           just = "top",
           x = unit(0.5, "lines"), 
           y = unit(7, "lines")
           ),
  textGrob(label = "B", 
           just = "top",
           x = unit(0.5, "lines"), 
           y = unit(18, "lines"))
  ),
  layout_matrix = rbind(
    c(5, 1, 3),
    c(6, 2, 4)),
  heights = c(2, 5),
  widths = c(.2, 5, 2))

# p_combined_color
ggsave(png_map_teo_color, 
       plot = p_combined_color, 
       device = "png", 
       width = 7.2, 
       height = 5,
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
                        theme(legend.position = "left") +
                        labs(fill = "Teosinte ranges")),
  cowplot::get_legend(p_map),
  textGrob(label = "A", 
         just = "top",
         x = unit(0.5, "lines"), 
         y = unit(7, "lines")
  ),
  textGrob(label = "B", 
         just = "top",
         x = unit(0.5, "lines"), 
         y = unit(18, "lines"))
  ),
  layout_matrix = rbind(
    c(5, 1, 3),
    c(6, 2, 4)),
  heights = c(2, 5),
  widths = c(.2, 5, 2))

#p_combined_black
ggsave(png_map_teo_black, 
       plot = p_combined_black, 
       device = "png", 
       width = 7.2,
       height = 5,
       units = "in",
       dpi = 300)
