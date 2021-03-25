#!/usr/bin/env Rscript
# working directory is hilo/
# calculation/plotting global ancestry for low-coverage maize/mexicana individuals 
# in HILO adaptation project
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

# load variables from snakemake
# get output plot and table filenames
png_structure = snakemake@output[["png_structure"]]
# png_structure = "global_ancestry/plots/HILO_MAIZE55_PARV50_structure_K3.png"
txt_allo_pops_alphas = snakemake[["allo_pops_alphas"]]
# txt_allo_pops_alphas = "global_ancestry/results/NGSAdmix/HILO_MAIZE55_PARV50/K3_alphas_by_allo_pop.txt"
txt_groups_alphas = snakemake[["groups_alphas"]]
# txt_groups_alphas = "global_ancestry/results/NGSAdmix/HILO_MAIZE55_PARV50/K3_alphas_by_group.txt"


# get colors for plot
# source("colors.R")
source(snakemake@input[["colors"]])

# get NGSAdmix results (d_admix2)
load(snakemake@input[["k3"]])
# load("global_ancestry/results/NGSAdmix/HILO_MAIZE55_PARV50/K3_alphas_by_ind.RData")

# summarise mean ancestry per group
options(scipen=999) # don't use scientific notation when printing alphas
d_admix2 %>%
  dplyr::select(group, LOCALITY, ELEVATION, maize, parviglumis, mexicana) %>%
  group_by(., group) %>%
  summarise(
    maize = mean(maize),
    parviglumis = mean(parviglumis),
    mexicana = mean(mexicana),
    n = n(),
    .groups = "drop") %>%
  write.table(.,
            txt_groups_alphas,
            col.names = T, row.names = F, quote = F, sep = "\t")
d_admix2 %>%
  dplyr::select(group, LOCALITY, ELEVATION, maize, parviglumis, mexicana) %>%
  dplyr::filter(group %in% c("allopatric_maize", "allopatric_mexicana", "parviglumis")) %>%
  group_by(., group, LOCALITY) %>%
  summarise(
    maize = mean(maize),
    parviglumis = mean(parviglumis),
    mexicana = mean(mexicana),
    n = n(),
    .groups = "drop") %>%
  write.table(.,
              txt_allo_pops_alphas,
              col.names = T, row.names = F, quote = F, sep = "\t")
options(scipen=0)

# make STRUCTURE-like ancestry plots:
p_structure_symp <- d_admix2 %>%
  filter(symp_allo == "sympatric") %>%
  arrange(., ELEVATION, popN) %>%
  group_by(zea, LOCALITY) %>% 
  mutate(sample = row_number()) %>%
  ungroup() %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize", "parviglumis")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_grid(zea ~ reorder(paste0(LOCALITY, "\n", "(", ELEVATION, "m", ")"), 
                                  ELEVATION), 
             scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_maize_mex_parv) +
  labs(fill = "Ancestry", x = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(angle = 90)) +
  ylab("Admixture proportion")

# p_structure_symp

# allopatric mexicana/maize/parviglumis structure only
p_structure_allo <- d_admix2 %>%
  filter(group %in% c("allopatric_maize", "allopatric_mexicana", "parviglumis")) %>%
  mutate(ELEVATION = ifelse(group == "allopatric_maize", 983, ELEVATION)) %>%
  mutate(ELEVATION = ifelse(group == "parviglumis", 1008, ELEVATION)) %>%
  arrange(., ELEVATION, popN) %>%
  mutate(sample = 1:nrow(.)) %>%
  mutate(sample = ifelse(group == "parviglumis", sample + 1*2,
                         ifelse(LOCALITY == "Puerta Encantada", sample + 2*2,
                                ifelse(LOCALITY == "Malinalco", sample + 3*2, 
                                       ifelse(LOCALITY == "Amecameca", sample + 4*2, 
                                              sample))))) %>%
  tidyr::gather(., "ancestry", "p", c("mexicana", "maize", "parviglumis")) %>%
  ggplot(., aes(fill = ancestry, y = p, x = sample)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = col_maize_mex_parv, labels = c("maize ancestry", "mexicana ancestry", "parviglumis ancestry")) +
  labs(fill = NULL, 
       y = "Admixture Proportion") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(name = element_blank(), 
                   limits = c(28, 25.5 + 57, 5.5 + 57 + 52, 21 + 57 + 52, 40 + 57 + 52),
                   #limits = c(27.5, 25 + 56, 5.5 + 56 + 51, 20 + 56 + 51, 37.5 + 56 + 51),
                   labels = c("Maize\nPalmar Chico\n(983m)",
                              "Parviglumis\nPalmar Chico\n(1008m)",
                              "Mexicana\nP. Encantada\n(1658m)", 
                              "Mexicana\nMalinalco\n(1887m)", 
                              "Mexicana\nAmecameca\n(2467m)")) +
  theme(plot.margin = margin(c(0,0,0,5.5)))
# p_structure_allo

# p_structure_allo <- d_admix2 %>%
#   filter(group %in% c("allopatric_maize", "allopatric_mexicana", "parviglumis")) %>%
#   mutate(ELEVATION = ifelse(LOCALITY == "Palmar Chico", 1008, ELEVATION)) %>%
#   arrange(., ELEVATION, popN) %>%
#   group_by(zea, LOCALITY) %>% 
#   mutate(sample = row_number()) %>%
#   ungroup() %>%
#   tidyr::gather(., "ancestry", "p", c("mexicana", "maize", "parviglumis")) %>%
#   ggplot(., aes(fill = ancestry, y = p, x = sample)) +
#   geom_bar(stat = "identity", position = "fill") + 
#   facet_grid(zea ~ reorder(LOCALITY, ELEVATION), scales = "free_x", space = "free_x") +
#   scale_fill_manual(values = col_maize_mex_parv) +
#   labs(fill = "Ancestry", x = "Sample") +
#   theme_classic() +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_text(angle = 90)) +
#   ylab("Admixture proportion")
# p_structure_allo

# multipanel plot of global ancestry data:
p_structure <- grid.arrange(grobs = list(ggplotGrob(p_structure_allo +
                                                     theme(
                                                       axis.text.x = element_text(size = 8),
                                                       legend.position = "None")),
                                        ggplotGrob(p_structure_symp +
                                                     theme(legend.position = "None")),
                                        cowplot::get_legend(p_structure_allo +
                                                              theme(legend.position = "bottom")),
                                        textGrob(label = "A", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0, "lines")),
                                        textGrob(label = "B", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0, "lines"))
                                       ),
                                        layout_matrix = rbind(
                                          c(4),
                                          c(3),
                                          c(1),
                                          c(5),
                                          c(2)),
                           heights = c(.1, .6, 3, .1, 7),
                           widths = c(1))

# p_structure
ggsave(png_structure, 
       plot = p_structure, 
       device = "png", 
       width = 7.5, 
       height = 5.5, 
       units = "in",
       dpi = 300)
