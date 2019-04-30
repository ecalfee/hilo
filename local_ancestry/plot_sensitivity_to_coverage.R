library(dplyr)
library(ggplot2)

# in this script I plot sensitivity to depth of coverage
# by taking one individual per population with >1.5x coverage
# and downsampling their read counts to lower coverage

depths <- c(0.05, 0.1, 0.25, 0.5, 1, 1.5)
IDs <- c("HILO207", "HILO41", "HILO218", "HILO42") # samples I downsampled
pops <- c("pop21", "pop26", "pop371", "pop370")

# what was the original coverage?
# what is the current coverage of each sample in IDs_to_sample?
original_file <- "../global_ancestry/results/NGSAdmix/pass2_alloMAIZE/globalAdmixtureByIncludedIndividual.txt"
original_coverage <- read.table(original_file, stringsAsFactors = F, header = T) %>%
  left_join(data.frame(ID = IDs, stringsAsFactors = F), ., by = "ID")


# GRAHAM's suggestion was also to mask sites per individual where they do not have any high-confidence call (maybe low information there)
# so, for example, I could look at just the sites in the genome where one of their posterior's is > .85
post0 <- read.table(paste0("results/ancestry_hmm/pass2_alloMAIZE/output_noBoot/", IDs[2], ".posterior"),
                    stringsAsFactors = F, header = T) # original
post1 <- read.table(paste0("results/ancestry_hmm/downsampled_1x/output_noBoot/", IDs[2], ".posterior"),
                    stringsAsFactors = F, header = T)

cor(post0$X2.0, post1$X2.0)
summary(lm(post0$X2.0 ~ post1$X2.0))

# summarize how well ancestry output of hmm with full coverage compares to ancestry output of hmm with downsampled coverage
results <- lapply(1:4, function(i){
  pop_ids <- read.table(paste0("results/ancestry_hmm/pass2_alloMAIZE/input/", pops[i], ".anc_hmm.ids"))$V1
  
  anc0 <- read.table(paste0("results/ancestry_hmm/pass2_alloMAIZE/output_noBoot/anc/", pops[i], ".anc.ind"),
                     stringsAsFactors = F, header = F)[ , which(pop_ids == IDs[i])] # original

  anc <- lapply(depths, function(d) 
    read.table(paste0("results/ancestry_hmm/downsampled_", d, "x/output_noBoot/anc/", pops[i], ".anc.ind"),
               stringsAsFactors = F, header = F)[ , which(pop_ids == IDs[i])])
  # plot difference between original and downsampled ancestry frequencies
  png(filename = paste0("plots/hist_diff_full_vs_downsampled_coverage_", IDs[i], ".png"),
      width = 12, height = 8, units = "in", res = 200)
  par(mfrow=c(2,3))
  sapply(1:6, function(j) hist(x = anc[[j]] - anc0, 
                               main = paste0("full x vs.", depths[j], "x ", IDs[i])))
  dev.off()
  
  return(list(predict_anc = lapply(1:6, function(j) summary(lm(anc0 ~ anc[[j]]))), # switched this to summary not full lm()
              summary = data.frame(corr_anc = sapply(1:6, function(j) cor(anc0, anc[[j]])),
              mean_anc = sapply(1:6, function(j) mean(anc[[j]])),
              depth = depths,
              pop = pops[i],
              ID = IDs[i], stringsAsFactors = F)))
})

summary <- do.call(bind_rows, 
             lapply(results, function(x) x[["summary"]]))

summary %>%
  ggplot(aes(x = depth, y = mean_anc, color = paste(ID, pop, sep = "_"))) +
  labs(color = "downsampled individual") +
  geom_line() +
  geom_point() +
  ggtitle("mean mexicana ancestry across all local ancestry calls")
ggsave("plots/alpha_est_ancestry_hmm_fullx_vs_downsampled.png", 
       height = 5, width = 8, units = "in", device = "png")

summary %>%
  ggplot(aes(x = depth, y = corr_anc, color = paste(ID, pop, sep = "_"))) +
  labs(color = "downsampled individual") +
  geom_line() +
  geom_point() +
  ggtitle("correlation in est. ancestry freq. between full x and downsampled x")
ggsave("plots/correlation_ancestry_hmm_fullx_vs_downsampled.png", 
       height = 5, width = 8, units = "in", device = "png")

summary %>%
  mutate(., r_squared = unlist(do.call(c, 
                                       lapply(results, function(x) 
                                         lapply(x$predict_anc, function(b) 
                                           summary(b)$r.squared))))) %>%
  ggplot(., aes(x = depth, y = r_squared, color = paste(ID, pop, sep = "_"))) +
  labs(color = "downsampled individual") +
  geom_line() +
  geom_point() +
  ggtitle("R-squared lm ancestry full x ~ ancestry downsampled x")
ggsave("plots/r-squared_ancestry_hmm_fullx_vs_downsampled.png", 
       height = 5, width = 8, units = "in", device = "png")
  
  
  
  