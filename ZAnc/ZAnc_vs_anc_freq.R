library(scales)
d <- read.table("results/pass1_maize/pop.ZAnc", stringsAsFactors = F,
                header = T)
summary(d)
quantile(d$mean_anc, c(.5, .95, .99, .999))
quantile(d$ZAnc, c(.5, .95, .99, .999))
png(file = "plots/corr_ZAnc_and_mean_ancestry.png",
    width = 10, height = 10, units = "in", res = 300)
with(d, plot(mean_anc, ZAnc, 
             cex = 1,
             col = alpha(ifelse(mean_anc > 0.55 & ZAnc > 13, "red", "blue"), .005),
             pch = 18,
             main = "99% ZAnc and Mean Ancestry across pops"))
abline(h = mean(d$ZAnc), v = mean(d$mean_anc))
dev.off()

mean(d$mean_anc)
sd(d$mean_anc)

# to do: what does the ancestry cov. matrix look like?
# for weird low ZAnc but moderate mean_anc (0.17 ish)
# plot what individual populations are doing

