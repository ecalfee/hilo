#!/usr/bin/env Rscript
# working directory is hilo/

library(dplyr)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)

Ne = snakemake@params[["Ne"]]
# Ne = 10000
yesno = snakemake@params[["yesno"]]
# yesno = "yes"

png_qq = snakemake@output[["png"]]
#png_qq = "ZAnc/plots/Ne10000_yesBoot/QQ.png"

# mean ancestry
# make dataframe with quantiles to plot
mex_maize = c("mexicana", "maize")
QQ_meanAnc = list(maize = NULL, mexicana = NULL)
for (zea in mex_maize){
  sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne", Ne, "_", yesno, "Boot/", zea, ".MVN.RData")
  anc_file = paste0("local_ancestry/results/ancestry_hmm/HILO_MAIZE55/Ne", Ne, "_", yesno, "Boot/anc/", zea, ".pops.anc.RData")
  load(anc_file)
  load(sim_file)
  QQ_meanAnc[[zea]] <- data.frame(p = seq(0, 1, length.out = 10000)) %>%
    mutate(q_MVN = quantile(MVN_mean, p), 
           q_data = quantile(anc_mean, p),
           zea = zea,
           group = paste("sympatric", zea))
}

# QQ plots:
p_meanAnc <- ggplot(data = do.call(bind_rows, QQ_meanAnc)) + 
    geom_point(aes(x = q_MVN, y = q_data, color = p)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Simulated mexicana ancestry") +
    ylab("Observed mexicana ancestry") +
    scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
    theme_light() +
    coord_fixed() +
    xlim(0:1) +
    ylim(0:1) +
    facet_wrap(~group)

# slope ancestry ~ elevation
# make dataframe with quantiles to plot
QQ_lmElev = list(maize = NULL, mexicana = NULL)
for (zea in mex_maize){
  sim_file = paste0("ZAnc/results/HILO_MAIZE55/Ne", Ne, "_", yesno, "Boot/", zea, ".lmElev.sim.RData")
  fit_file = paste0("ZAnc/results/HILO_MAIZE55/Ne", Ne, "_", yesno, "Boot/", zea, ".lmElev.fit.RData")
  load(fit_file)
  load(sim_file)
  QQ_lmElev[[zea]] <- data.frame(p = seq(0, 1, length.out = 10000)) %>%
    mutate(q_MVN = quantile(fits_sim$envWeights, p), 
           q_data = quantile(fits$envWeights, p),
           zea = zea,
           group = paste("sympatric", zea))
}

# QQ plots:
p_lmElev <- ggplot(data = do.call(bind_rows, QQ_lmElev)) + 
  geom_point(aes(x = q_MVN, y = q_data, color = p)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Simulated slope (mexicana anc ~ elev)") +
  ylab("Observed slope") +
  scale_color_viridis(option = "plasma", direction = -1, name = "Quantile") +
  theme_light() +
  coord_fixed() +
  scale_y_continuous(limits = with(do.call(bind_rows, QQ_lmElev), c(min(c(q_MVN, q_data)), max(c(q_MVN, q_data)))),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(limits = with(do.call(bind_rows, QQ_lmElev), c(min(c(q_MVN, q_data)), max(c(q_MVN, q_data)))),
                     labels = scales::number_format(accuracy = 0.01)) +
  facet_wrap(~group)

p_combined <- grid.arrange(grobs = list(textGrob(label = "A", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0.5, "lines")),
                                        ggplotGrob(p_meanAnc),
                                        textGrob(label = "B", 
                                                 x = unit(0.5, "lines"), 
                                                 y = unit(0.5, "lines")),
                                        ggplotGrob(p_lmElev)),
                           nrow = 4,
                           heights = c(.2, 2, .1, 2)
)
#p_combined
ggsave(png_qq, 
       plot = p_combined, 
       device = "png", 
       width = 6, 
       height = 6, 
       units = "in",
       dpi = 300)
