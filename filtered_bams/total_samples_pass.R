#!/usr/bin/env Rscript

# this script summarises filtered bam metrics and creats a
# plot of samples that pass minimum coverage threshold
# and lists of samples to use as input for other analyses
library(dplyr)
library(tidyr)
library(ggplot2)

#colors_file = snakemake@input[["colors"]]
#meta_file = snakemake@input[["colors"]]
#colors_file = "../colors.R" 
#meta_file = "../samples/HILO_MAIZE55_meta.RData"
setwd(snakemake@params[["wd"]]) # set working directory to the path where this script is
source("../colors.R") # colors for plots
load("../samples/HILO_MAIZE55_meta.RData") # sample metadata, including estimated sequence coverage

# plot total included samples per population
p_seq_counts <- meta %>%
  filter(dataset == "HILO") %>%
  mutate(coverage = cut(est_coverage, 
                        breaks = c(0, 0.05, 0.5, 10)
  )) %>%
  ggplot(., aes(alpha = coverage,
                fill = zea,
                x = LOCALITY)) +
  geom_bar() +
  labs(fill = "Zea", x = "Population", y = "# Individuals") +
  scale_alpha_discrete(range = c(.6, 1),
                       name = "Coverage",
                       labels = c(expression("x "<" 0.5"), expression("x ">=" 0.5"))) +
  scale_fill_manual(values = col_maize_mex_parv) +
  facet_grid(zea~.) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))#
#p_seq_counts
ggsave("../../hilo_manuscript/figures/p_seq_counts.png",
       plot = p_seq_counts,
       height = 4, width = 5.4, units = "in",
       device = "png")
ggsave("plots/p_seq_counts.png",
       plot = p_seq_counts,
       height = 4, width = 5.4, units = "in",
       device = "png")

# write ID and bam files by group/population
# all individuals included
meta_bams <- meta %>%
  mutate(bam = ifelse(dataset == "HILO",
                      paste0("filtered_bams/merged_bams/", ID, ".sort.dedup.bam"),
                      paste0("filtered_bams/Maize55/", ID, ".sort.dedup.bam")))
# allopatric_maize, sympatric_maize etc.
for (g in unique(meta_bams$group)){
  filter(meta_bams, group == g) %>%
    dplyr::select(ID) %>%
    write.table(., paste0("../samples/ALL_byPop/", g, "_ids.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
  filter(meta_bams, group == g) %>%
    dplyr::select(bam) %>%
    write.table(., paste0("../samples/ALL_byPop/", g, "_bams.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
}
# pop22, pop366 etc.
for (p in unique(meta_bams$popN)){
  filter(meta_bams, popN == p) %>%
    dplyr::select(ID) %>%
    write.table(., paste0("../samples/ALL_byPop/pop", p, "_ids.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
  filter(meta_bams, popN == p) %>%
    dplyr::select(bam) %>%
    write.table(., paste0("../samples/ALL_byPop/pop", p, "_bams.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
}

# all individuals included with over 0.5x coverage
# (i.e. used in local ancestry inference)
# allopatric_maize, sympatric_maize etc.
for (g in unique(meta_bams$group)){
  filter(meta_bams, group == g, est_coverage >= 0.5) %>%
    dplyr::select(ID) %>%
    write.table(., paste0("../samples/Over0.5x_byPop/", g, "_ids.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
  filter(meta_bams, group == g, est_coverage >= 0.5) %>%
    dplyr::select(bam) %>%
    write.table(., paste0("../samples/Over0.5x_byPop/", g, "_bams.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
}
# pop22, pop366 etc.
for (p in unique(meta_bams$popN)){
  filter(meta_bams, popN == p, est_coverage >= 0.5) %>%
    dplyr::select(ID) %>%
    write.table(., paste0("../samples/Over0.5x_byPop/pop", p, "_ids.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
  filter(meta_bams, popN == p, est_coverage >= 0.5) %>%
    dplyr::select(bam) %>%
    write.table(., paste0("../samples/Over0.5x_byPop/pop", p, "_bams.list"), 
                col.names = F, row.names = F,
                sep = "\t", quote = F)
}

# ------------------------------------------
# what is a reasonable cut-off for % of individuals with coverage to call a SNP?
# under a basic poisson model, what % of individuals typically have data for any SNP?
# just HILO
included <- filter(meta, dataset == "HILO")
sims <- sapply(1:10000, function(i) sum(rpois(n = length(included$est_coverage), 
                                              lambda = included$est_coverage) > 0)/length(included$est_coverage))
# with an estimated loss of coverage for not meeting quality > 20
sims2 <- sapply(1:10000, function(i) sum(rpois(n = length(included$est_coverage), 
                                               lambda = included$est_coverage*2/3) > 0)/length(included$est_coverage))
# with overdispersion (negative binomial)
sims3 <- sapply(1:10000, function(i) sum(rnbinom(n = length(included$est_coverage),
                                                 size = 1,
                                                 mu = included$est_coverage*2/3) > 0)/length(included$est_coverage))

# also including MAIZE
sims_w55 <- (sims*length(included$est_coverage)+55)/(length(included$est_coverage) + 55)

sims2_w55 <- (sims2*length(included$est_coverage)+55)/(length(included$est_coverage) + 55)

sims3_w55 <- (sims3*length(included$est_coverage)+55)/(length(included$est_coverage) + 55)
# table(sims2_w55 > 0.5)/length(sims2_w55)
# about 3.5% of the genome approximately has < 50% coverage across sampled individuals
# if we assume poisson. But with overdispersion it's potentially a lot more
# hist(sims3_w55)
# > 40% of individuals should be most of the genome that has regular mapping even with overdispersion