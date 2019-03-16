library(dplyr)
library(ggplot)
library(reader)

getTheta <- function(theta_file, pop, region){
  if (file.exists(theta_file) && reader::file.nrow(theta_file) > 1){
    d <- read.table(theta_file, sep = "\t")[ , c(2:5,14)]
    colnames(d) <- c("chr", "pos", "theta_W_sum", "theta_pi_sum", "n_bp")
    d1 <- d %>% # scale thetas per bp not per region analyzed, based on total number of bp's with data
        mutate(theta_W = theta_W_sum/n_bp) %>%
        mutate(theta_pi = theta_pi_sum/n_bp) %>%
        mutate(pop = pop) %>%
        mutate(region = region)
    
  }else{ # returns empty dataframe if no file exists (some regions some populations have no data)
    d1 <- NULL
    print(paste("warning: no theta file found for pop", pop, "region", region))
  }
  return(d1)
}

getFst <- function(fst_file, pop, region){
  if (file.exists(theta_file)){
  d <- read.table(fst_file, sep = "\t")
  colnames(d) <- c("fst_Hudson", "fst_avg_ratios")
  d1 <- d %>%
    mutate(pop = pop) %>%
    mutate(region = region)
  }else{
    d1 <- NULL
    print(paste("warning: no fst file found for pop", pop, "region", region))
  }
  return(d1)
}

#name = "outliers"
#name = "outliers_100kb_flanking"
#name = "outliers_mean_anc0.6"
name = "outliers_mean_anc0.6_100kb_flanking"

outliers <- read.table(paste0("../ZAnc/results/pass1_maize/", name, ".bed"), stringsAsFactors = F, header = F)

colnames(outliers) <- c("chr", "start", "end", "max_ZAnc", "max_pop_freq")
fst_pairs <- c("symp.maize-allo.mexicana", "symp.maize-symp.mexicana", "symp.mexicana-allo.mexicana")

pops <- unique(read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                          header = T, sep = "\t")[ , c("popN", "zea", "symp_allo", "LOCALITY")]) %>%
  mutate(group = paste(substr(symp_allo, 1, 4), zea, sep = ".")) %>%
  mutate(pop = paste0("pop", popN))
groups <- unique(pops$group)
unique(pops$popN)
pops_all <- sort(pops$pop)
#pops_incl <- pops_all[c(1:15, 17, 22:25, 27, 29:30)] # some populations don't have proper files (maybe no homozygotes?)

regions <- with(outliers, paste0("chr", chr, "_", start, "-", end))
regions_incl <- regions[-c(1)] # missing Fst stats for this region -- skip for now

thetas <- do.call(rbind,
                  lapply(regions, function(r) do.call(rbind,
                  lapply(groups, function(pop) getTheta(paste0("results/piFst_", name, "/mex2/", r, "/", pop, ".thetasAll.pestPG"), 
                         pop = pop, region = r)))))


ggplot(thetas, aes(x = pop, y = theta_pi, color = region)) +
  geom_point() + 
  geom_line(aes(group = region)) +
  ggtitle(paste0("pi across groups for high mex ", name))
ggsave(filename = paste0("plots/pi_groups_", name, ".png"), device = "png",
       plot = last_plot(), 
       width = 10, height = 10, units = "in",
       dpi = 300)


FSTs <- do.call(rbind,
                #lapply(regions_incl, function(r) do.call(rbind,
                lapply(regions, function(r) do.call(rbind,
                                                      lapply(fst_pairs, function(pop) 
                                                        getFst(paste0("results/piFst_", name, "/mex2/", r,
                                                                         "/", pop, ".fst.stats"), pop = pop, region = r)))))

ggplot(FSTs, aes(x = pop, y = fst_Hudson, color = region)) +
  geom_point() + 
  geom_line(aes(group = region)) +
  ggtitle(paste0("fst across groups for high mex ", name))
ggsave(filename = paste0("plots/fst_groups_", name, ".png"), device = "png",
       plot = last_plot(), 
       width = 10, height = 10, units = "in",
       dpi = 300)



thetasWind <- do.call(rbind,
                  lapply(regions, function(r) do.call(rbind, 
                                                      lapply(groups, function(pop) 
                                                        getTheta(paste0("results/piFst_", name, "/mex2/", r, "/", pop, ".thetasWindows.pestPG"), 
                                                                 pop = pop, region = r)))))
                        
for (j in regions){
  p <- thetasWind %>%
    filter(region==j) %>%
    ggplot(., aes(x = pos, y = theta_pi, color = pop)) +
    geom_point() + 
    ggtitle(paste0("pi for high mex ", name, ": ", j))
  ggsave(filename = paste0("plots/pi_", name, "_region_", j, ".png"), 
         device = "png",
         plot = p, 
         width = 10, height = 5, units = "in",
         dpi = 300)
}


thetas_ind <- do.call(rbind,
                  lapply(regions, function(r) do.call(rbind,
                                                      lapply(unique(pops$pop), function(pop) getTheta(paste0("results/piFst_", name, "/mex2/", r, "/", pop, ".thetasAll.pestPG"), 
                                                                                          pop = pop, region = r)))))
# compare pi across groups to pi across pops for outlier regions:
left_join(pops, thetas_ind) %>%
  ggplot(., aes(x = region, y = theta_pi, color = group)) +
  geom_point() + 
  geom_line(aes(group = pop)) +
  ggtitle(paste0("pi across pops for high mex ", name))
ggsave(filename = paste0("plots/pi_ind_pops_", name, ".png"), 
       device = "png",
       width = 10, height = 5, units = "in",
       dpi = 300)
left_join(pops, thetas_ind) %>%
  filter(group == "symp.maize") %>%
  ggplot(., aes(x = region, y = theta_pi, color = LOCALITY)) +
  geom_point() + 
  geom_line(aes(group = pop)) +
  ggtitle(paste0("pi across symp. maize pops for high mex ", name))
left_join(pops, thetas_ind) %>%
  ggplot(., aes(x = region, y = theta_pi, color = LOCALITY)) +
  geom_point(aes(group = pop)) + 
  facet_wrap(~group) +
  ggtitle(paste0("pi across locality pops for high mex ", name))
ggsave(filename = paste0("plots/pi_locality_pops_", name, ".png"), 
       device = "png",
       width = 10, height = 5, units = "in",
       dpi = 300)
left_join(pops, thetas_ind) %>%
  group_by(group, region) %>%
  summarize(theta_pi = mean(theta_pi)) %>%
  ggplot(., aes(x = region, y = theta_pi, color = group)) +
  geom_point() + 
  geom_line(aes(group = group)) +
  ggtitle(paste0("pi mean across pops for high mex ", name))
ggsave(filename = paste0("plots/pi_mean_pops_", name, ".png"), 
       device = "png",
       width = 10, height = 5, units = "in",
       dpi = 300)
ggplot(thetas, aes(x = region, y = theta_pi, color = pop)) +
  geom_point() + 
  geom_line(aes(group = pop)) +
  ggtitle(paste0("pi across groups of ind. for high mex ", name))
ggsave(filename = paste0("plots/pi_grouped_inds_", name, ".png"), 
       device = "png",
       width = 10, height = 5, units = "in",
       dpi = 300)
#pops_incl[c(1:10, 12, 14:16, 18:22)]
thetasWind_ind <- do.call(rbind,
                      lapply(regions, function(r) do.call(rbind, 
                                                          lapply(unique(pops$pop), function(pop) 
                                                            getTheta(paste0("results/piFst_", name, "/mex2/", r, "/", pop, ".thetasWindows.pestPG"), 
                                                                     pop = pop, region = r)))))

#d <- read.table(paste0("results/piFst_", name, "/mex2/", regions[1], "/", unique(pops$pop)[5], ".thetasWindows.pestPG"), sep = "\t") 


for (j in unique(thetasWind_ind$region)){
  p <- left_join(thetasWind_ind, pops, by = "pop") %>%
    filter(region==j) %>%
    ggplot(., aes(x = pos/1000, y = theta_pi, color = LOCALITY)) +
    xlab("pos (kb)") +
    geom_line() + 
    ggtitle(paste0("pi for high mex ", name, ": ", j)) +
    facet_wrap(~group)
  plot(p)
  ggsave(filename = paste0("plots/pi_ind_pop_", name, "_region_", j, ".png"), 
         device = "png",
         plot = p, 
         width = 10, height = 5, units = "in",
         dpi = 300)
}


