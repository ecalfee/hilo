library(dplyr)
library(ggplot)

getTheta <- function(theta_file, pop, region){
  d <- read.table(theta_file)[ , c(2:5,14)]
  colnames(d) <- c("chr", "pos", "theta_W_sum", "theta_pi_sum", "n_bp")
  d1 <- d %>% # scale thetas per bp not per region analyzed, based on total number of bp's with data
    mutate(theta_W = theta_W_sum/n_bp) %>%
    mutate(theta_pi = theta_pi_sum/n_bp) %>%
    mutate(pop = pop) %>%
    mutate(region = region)
  return(d1)
}
getFst <- function(fst_file, pop, region){
  d <- read.table(fst_file)
  colnames(d) <- c("fst_Hudson", "fst_avg_ratios")
  d1 <- d %>%
    mutate(pop = pop) %>%
    mutate(region = region)
  return(d1)
}



exclude_slurm_tasks <- c(0,6) # errors running
#outliers <- read.table("../ZAnc/results/pass1_maize/outliers.bed", stringsAsFactors = F, header = F)[-(exclude_slurm_tasks+1), ]

#colnames(outliers) <- c("chr", "start", "end", "max_ZAnc", "max_pop_freq")
fst_pairs <- c("symp.maize-allo.mexicana", "symp.maize-symp.mexicana", "symp.mexicana-allo.mexicana")
pops <- c("symp.maize", "symp.mexicana", "allo.mexicana") 
#i = 1
#outliers[i,]
#dir <- paste0("results/piFst_outliers/mex2/chr", outliers[i, "chr"], "_", outliers[i,"start"], "-", outliers[i, "end"])
regions <- read.table("results/piFst_outliers/mex2/regions.txt", stringsAsFactors = F, header = F)[-(exclude_slurm_tasks+1), ]

thetas <- do.call(rbind,
                  lapply(regions, function(r) do.call(rbind,
                  lapply(pops, function(pop) getTheta(file.path("results/piFst_outliers/mex2", r,
                                                                paste0(pop, ".thetasAll.pestPG")), pop, region = r)))))



#thetas <- do.call(rbind,
#                  lapply(1:nrow(outliers), function(r) do.call(rbind,
#                                                               lapply(pops, function(pop) 
#                                                                 file.exists(file.path(paste0("results/piFst_outliers/mex2/chr", outliers[r, "chr"], "_", outliers[i,"start"], "-", outliers[i, "end"]),
#                                                                                                             paste0(pop, ".thetasAll.pestPG")))))))





ggplot(thetas, aes(x = pop, y = theta_pi, color = region)) +
  geom_point() + 
  geom_line(aes(group = region)) +
  ggtitle("pi across groups for high mexicana ancestry ZAnc outlier regions")
ggsave(filename = "pi_outliers.png", device = "png",
       plot = last_plot(), 
       width = 10, height = 10, units = "in",
       dpi = 300)


FSTs <- do.call(rbind, # some excluded (file issues)
                  lapply(regions[c(1:3,6:26)], function(r) do.call(rbind,
                                                      lapply(fst_pairs, function(pop) getFst(file.path("results/piFst_outliers/mex2", r,
                                                                                                       paste0(pop, ".fst.stats")), pop, region = r)))))
ggplot(FSTs, aes(x = pop, y = fst_Hudson, color = region)) +
  geom_point() + 
  geom_line(aes(group = region)) +
  ggtitle("fst across groups for high mexicana ancestry ZAnc outlier regions")
ggsave(filename = "fst_outliers.png", device = "png",
       plot = last_plot(), 
       width = 10, height = 10, units = "in",
       dpi = 300)



thetasWind <- do.call(rbind,
                  lapply(regions[c(1:3,6:26)], function(r) do.call(rbind,
                                                      lapply(pops, function(pop) getTheta(file.path("results/piFst_outliers/mex2", r,
                                                                                                    paste0(pop, ".thetasWindows.pestPG")), pop, region = r)))))
for (j in regions){
  p <- thetasWind %>%
    filter(region==j) %>%
    ggplot(., aes(x = pos, y = theta_pi, color = pop)) +
    geom_point() + 
    ggtitle(paste0("pi diversity for high mexicana ancestry ZAnc outlier: ", j))
  ggsave(filename = paste0("pi_region_", j, ".png"), device = "png",
         plot = p, 
         width = 10, height = 5, units = "in",
         dpi = 300)
}




do.call(rbind,
                lapply(regions, function(r) do.call(rbind,
                                                    lapply(fst_pairs, function(pop) file.exists(file.path("results/piFst_outliers/mex2", r,
                                                                                                     paste0(pop, ".fst.stats")))))))
do.call(rbind,
        lapply(regions, function(r) do.call(rbind,
                                            lapply(fst_pairs, function(pop) with(file.info(file.path("results/piFst_outliers/mex2", r,
                                                                                                  paste0(pop, ".fst.stats"))), size <18)))))


lapply(fst_pairs, function(pop) getFst(file.path("results/piFst_outliers/mex2", "chr7_168608960-169087801",
                 paste0(pop, ".fst.stats")), pop = pop, region ="chr7_168608960-169087801"))


thetasWind <- do.call(rbind,
                  lapply(1:nrow(outliers), function(r) do.call(rbind,
                                                               lapply(pops, function(pop) getTheta(file.path(paste0("results/piFst_outliers/mex2/chr", outliers[i, "chr"], "_", outliers[i,"start"], "-", outliers[i, "end"]),
                                                                                                             paste0(pop, ".thetasWindows.pestPG")), pop, region = r)))))



thetasWind <- do.call(rbind,
                  lapply(pops, function(pop) getTheta(file.path(dir, paste0(pop, ".thetasWindows.pestPG")), pop, region)))
                 