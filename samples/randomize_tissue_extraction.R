planted <- read.csv("seeds_planted_5.1.19.csv")
set.seed(200)
planted_randomized <- planted[sample(1:nrow(planted)), ]
failed <- read.table("2019_25_May_failed2grow.txt", stringsAsFactors = F, header = F, sep = " ")
colnames(failed) <- c("RIMMEA", "popN", "family")
planted_randomized$failed <- paste0(planted_randomized$popN, substr(planted_randomized$family, 1, 1)) %in% paste0(failed$popN, substr(failed$family, 1, 1))
get_tissue <- planted_randomized[!planted_randomized$failed, ]
get_tissue$number <- rep(1:12, 16)[1:nrow(get_tissue)] # add well numbers
get_tissue$letter <- rep(sapply(LETTERS[1:8], function(x) rep(x, 12)), 2)[1:nrow(get_tissue)]

write.table(get_tissue, "get_tissue_5.28.19.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")

