library(NMF)
data = read.table("NMF_excluding_peaks_in_1kbTSS/data/mat.txt")
dim(data)
res1 <- nmf(data, 2:15, nrun = 100, .options = "pv", seed = 123456)
# evaluate rank survey
pdf("NMFRankSurvey.pdf")
plot(res1)
dev.off()

# save an object to a file
saveRDS(res1, file = "res1.rds")
