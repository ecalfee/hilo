# this script takes in an ancestry_hmm input file
# IDs, and coverage for included individuals in that pop
# and dowsamples specified higher coverage individuals 
# to 1.5x, 1x, .5x, .25x, and .1x coverage
library(dplyr)

IDs_to_downsample <- c("HILO207", "HILO41", "HILO218", "HILO42")

# what is the current coverage of each sample in IDs_to_sample?
coverage_file <- "../global_ancestry/results/NGSAdmix/pass2_alloMAIZE/globalAdmixtureByIncludedIndividual.txt"
original_coverage <- read.table(coverage_file, stringsAsFactors = F, header = T) %>%
  left_join(data.frame(ID = IDs_to_downsample, stringsAsFactors = F), ., by = "ID")

# keep each observed read (major or minor allele) with probability
# p = target_coverage/original_coverage
downsample <- function(target, original, counts){
  p <- target/original
  new_counts <- sapply(counts, function(i)
    rbinom(n = 1, size = i, prob = p))
  return(new_counts)
}


dir_input <- "results/ancestry_hmm/pass2_alloMAIZE/input"
dir_output <- "results/ancestry_hmm"

# make output file for each population and each target coverage
for (popN in unique(original_coverage$popN)){
  # which individuals are included and part of this pop?
  ids_file <- paste0(dir_input, "/pop", popN, ".anc_hmm.ids")
  counts_file <- paste0(dir_input, "/pop", popN, ".anc_hmm.input")
  pop_ids <- read.table(ids_file, stringsAsFactors = F, header = F)$V1
  
  # which of these individuals are being downsampled?
  pop_ids_to_downsample <- pop_ids[pop_ids %in% IDs_to_downsample]
  
  # get original full data ancestry_hmm input file
  counts <- read.table(counts_file, stringsAsFactors = F, header = F)
  # add column names; samples are in same order as pop_ids
  colnames(counts) <- c("chr", "pos", "maize_A", "maize_a", "mex_A", "mex_a",
                        "rDist", paste(sapply(pop_ids, function(x) rep(x, 2)), sep = "_", rep(c("A", "a"), length(pop_ids))))
  
  for (target in c(1.5, 1, 0.5, 0.25, 0.1, 0.05)){
    counts2 <- counts
    # downsample counts
    for (id in pop_ids_to_downsample){ # downsample all individuals in pop
      for (a in c("A", "a")){ # downsample both alleles
        counts2[ , paste(id, a, sep = "_")] <- downsample(target = target, 
                                                                        original = original_coverage$est_coverage[original_coverage$ID == id],
                                                                        counts = counts[ , paste(id, a, sep = "_")])
      }
    }
    new_dir <- paste0(dir_output, "/downsampled_", target, "x/input")
    dir.create(new_dir)
    write.table(counts2, 
                paste0(new_dir, "/pop", popN, ".anc_hmm.input"),
                col.names = F, row.names = F, quote = F)
    }
}



