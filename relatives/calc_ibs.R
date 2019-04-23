# calculate identity by state (ibs) for an estimate of relatedness:
ibs_ind <- read.table("results/N1000.L100/all.ibs", stringsAsFactors = F, header = T)
ibs_pair <- read.table("results/N1000.L100/all.ibspair", stringsAsFactors = F, header = T)
colnames(ibs_pair)
all_geno <- c("AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT")
all_same_geno <- paste0("p", all_geno, "_", all_geno)
# what percent is inferred to be ibs 2?
ibs2 <- ibs_pair %>%
  mutate(ibs2 = apply(dplyr::select(ibs_pair, all_same_geno), 1, sum))
IDs <- data.frame(ID = do.call(c, lapply(read.table(paste0("../samples/", "duplicates", "_IDs.list"), header = F, stringsAsFactors = F)$V1,
             function(x) rep(x, 2))),
             ind = 0:73,
            stringsAsFactors = F)
# which have the same ID?
ibs3 <- dplyr::rename(IDs, ID1 = ID) %>%
  rename(., ind1 = ind) %>%
  left_join(ibs2, ., by = "ind1")


ibs4 <- dplyr::rename(IDs, ID2 = ID) %>%
  rename(., ind2 = ind) %>%
  left_join(ibs3, ., by = "ind2") %>%
  dplyr::mutate(same_ID = (ID1 == ID2)) %>%
  dplyr::mutate(., ibs_AA_GG_TT_CC = pAA_AA + pGG_GG + pTT_TT + pCC_CC) %>%
  dplyr::filter(., !is.na(.$ind1) | .$ind1=="NA")

ggplot(ibs4, aes(x = ID1, y = ibs2, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("2 alleles IBS")
ggsave("plots/duplicates_N1000.L100_ibs2_allgeno_only.png",
       device = "png",
       width = 5, height = 5, units = "in",
       dpi = 200)
ggplot(ibs4, aes(x = ID1, y = ibs_AA_GG_TT_CC, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("2 alleles IBS - homozygouts genotypes")
ggsave("plots/duplicates_N1000.L100_ibs_homozyg_only.png",
       device = "png",
       width = 5, height = 5, units = "in",
       dpi = 200)
ggplot(ibs4, aes(x = ibs2, y = ibs_AA_GG_TT_CC, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("correlation same genotype het vs. hom")
ggsave("plots/duplicates_N1000.L100_ibs2_homozyg_vs_allgeno.png",
       device = "png",
       width = 5, height = 5, units = "in",
       dpi = 200)
# why are there multiple perfect matches? low coverage??
ibs4[ibs4$ibs2==1 & complete.cases(ibs4), c("ind1", "ind2", "ID1", "ID2", "nSites", "Llike", "ibs2", "ibs_AA_GG_TT_CC")]
