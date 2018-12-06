# gets column numbers to pull out of full beagle.gz file
d = read.table("../data/HILO_IDs_cov_pass1.csv", header = T, 
               stringsAsFactors = F)
header = t(read.table("../data/pass1_bam.all.GLheader", 
                           stringsAsFactors = F))
inds_include = d[((d$zea == "mexicana" & d$symp_allo == "allopatric") | 
                            (d$zea == "maize" & d$symp_allo == "sympatric")) & 
                            d$est_coverage >= 0.05 & !is.na(d$est_coverage), "pass1_ID"]
cols_include = which(header %in% c("marker", "allele1", "allele2", inds_include))

write.table(paste(cols_include, collapse=","), "../data/pass1_bam.SympMaizeAlloMex05.GLcolN",
          quote = F, col.names = F, row.names = F)
