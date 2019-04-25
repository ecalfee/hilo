# calculate identity by state (ibs) for an estimate of relatedness:
#PREFIX = "N1000.L100"
#NAME = "all"
PREFIX = "N1000L10kb"
NAME = "duplicates"
ibs_ind <- read.table(paste0("results/", PREFIX, "/", NAME, ".ibs"), stringsAsFactors = F, header = T)
ibs_pair <- read.table(paste0("results/", PREFIX, "/", NAME, ".ibspair"), stringsAsFactors = F, header = T)
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

# I want to plot matches to both sets of, e.g. HILO80 (first set ID1, second set ID2):
ibs4 %>%
  gather(., "set", "ID", c("ID1", "ID2")) %>%
  ggplot(., aes(x = ID, y = ibs2, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("2 alleles IBS") +
  facet_wrap(~set) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0("plots/duplicates_", PREFIX, "_ibs2_allgeno_only.png"),
       device = "png",
       width = 15, height = 10, units = "in",
       dpi = 200)

# I want to plot matches to both sets of, e.g. HILO80 (first set ID1, second set ID2):
ibs4 %>%
  gather(., "set", "ID", c("ID1", "ID2")) %>%
  ggplot(., aes(x = ID, y = ibs_AA_GG_TT_CC, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("2 alleles IBS") +
  facet_wrap(~set) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("plots/duplicates_", PREFIX, "_ibs_homozyg_only.png"),
       device = "png",
       width = 15, height = 10, units = "in",
       dpi = 200)
ibs4 %>%
  gather(., "set", "ID", c("ID1", "ID2")) %>%
  ggplot(., aes(x = ibs2, y = ibs_AA_GG_TT_CC, color = same_ID, size = nSites)) +
  geom_point(alpha = .5) +
  ggtitle("correlation same genotype het vs. hom")
ggsave(paste0("plots/duplicates_", PREFIX, "_ibs2_homozyg_vs_allgeno.png"),
       device = "png",
       width = 5, height = 5, units = "in",
       dpi = 200)
# why are there multiple perfect matches? low coverage??
ibs4[ibs4$ibs2==1 & complete.cases(ibs4), c("ind1", "ind2", "ID1", "ID2", "nSites", "Llike", "ibs2", "ibs_AA_GG_TT_CC")]
ibs4 %>%
  group_by(., same_ID) %>%
  summarise(., mean_ibs2 = mean(ibs2, na.rm = T))
# only look at higher coverage pairs:
ibs4 %>% # makes little difference in means
  filter(nSites > 10000) %>%
  group_by(., same_ID) %>%
  summarise(., mean_ibs2 = mean(ibs2, na.rm = T))
