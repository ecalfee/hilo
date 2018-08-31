# this script plots posterior results from ancestry_hmm
library(dplyr)
library(tidyr)
library(ggplot2)

pass1 <- read.table("../data/pass1_ids.txt", stringsAsFactors = F, 
                                     header = T, sep = "\t")
dir_in = "../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/output/firstgo_pop366_job26291088/"
popN = 366
ids_list = read.table(paste0("../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input/pop",
                  popN, ".anc_hmm.ids"), stringsAsFactors = F)$V1
for (id in ids_list){
  coverage = round(pass1[pass1$n == id, "est_coverage"], 2)
  post = read.table(paste0(dir_in, "/HILO", id, ".posterior"), stringsAsFactors = F, header = T) %>%
    rename(., mex.mex = X0.2) %>%
    rename(., maize.mex = X1.1) %>%
    rename(., maize.maize = X2.0)
  post1 = tidyr::gather(post, "anc", "p", 3:5)
  small1 = tidyr::gather(post[c(T, rep(F, 99)), ], "anc", "p", 3:5)
  # percent sites inferred with high confidence
  # /3 because there are 3 rows (and hmm states) 
  # per SNP and only one 
  # can possibly be over 50% posterior
  over9 = round(sum(small1$p > .9)/(nrow(small1)/3),2)
  over8 = round(sum(small1$p > .8)/(nrow(small1)/3),2)
  # estimated global ancestry from hmm
  perc_mex_est = small1 %>% 
    group_by(., anc) %>%
    summarise(., alpha = mean(p)) %>%
    .$alpha %*% (0:2)/2
  perc_mex_est = round(perc_mex_est, 2)
  
  # plot mexicana ancestry genomewide
  p1 = small1 %>%
    filter(anc == "mex.mex") %>%
    ggplot(aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p1 + facet_wrap(~chrom)
  ggsave(filename = paste0("../plots/tractsMexAnc_HILO", id, "_allCHR"),
         plot = p1 + facet_wrap(~chrom), 
         device = "png", height = 12, width = 20, units = "in")
  
  # plot posterior for mexicana 1 or 2 alleles genomewide
  p2 = small1 %>%
    filter(anc != "maize.maize") %>%
    ggplot(aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p2 + facet_wrap(~chrom)
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_allCHR"),
         plot = p2 + facet_wrap(~chrom), 
         device = "png", height = 12, width = 20, units = "in")
  
  p3 = post1 %>%
    filter(., chrom == 3,
           position < 1000000) %>%
    ggplot(., aes(x=position, y = p)) +
    geom_point(aes(color = anc)) + 
    ggtitle(paste0("zoomed in posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p3
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_zoomCHR3"),
         plot = p3, 
         device = "png", height = 12, width = 20, units = "in")
    
  p4 = small1 %>%
    filter(., chrom == 4 &
           position > 1e8*1.6 &
             position < 1e8*1.9) %>%
    ggplot(., aes(x=position, y = p)) +
    geom_point(aes(color = anc)) +
    ggtitle(paste0("zoomed in posterior prob. HILO", id, " est. cov = ", coverage), 
            paste0("est: % mex=", perc_mex_est, " % sites > .9 post= ", over9))
  p4
  ggsave(filename = paste0("../plots/tractsMexMaizeAnc_HILO", id, "_zoomCHR4"),
         plot = p4, 
         device = "png", height = 12, width = 20, units = "in")
  
}


