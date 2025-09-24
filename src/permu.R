library(plotrix)

load("diffbind/results/diffbind_work_space.RData")
# dbObj$called (94157, 22): binary
# dbObj$binding (94157, 25): similar with dbObj$called but entires are normalied scores

# mat's nrow as number of permutation (of order of 22 samples)
mat = matrix(0, 1000, 22)
dat = dbObj$called
print(dim(dat))
for(r in 1:nrow(mat)){
    for(c in 1:ncol(mat)){
        if(c==1){
            mat[r,1] = sum(dat[,1])
        }
        else{
            mat[r, c] = sum(apply(dat[,1:c], 1 , any))
        }
    }
    dat = dat[, sample(ncol(dat))]
}
write.table(mat, "mat.txt", row.names = F)

# plot 95% confidence interval
alpha = 0.05
data = read.table("mat.txt", header = TRUE)
n = dim(data)[1]
dat = data.frame(x = 1:dim(data)[2],
                 y = apply(data, 2, mean),
                 low = apply(data, 2, mean) - qt(1 - alpha/2, n - 1) * apply(data, 2, sd)/sqrt(n),
                 up = apply(data, 2, mean) + qt(1 - alpha/2, n - 1) * apply(data, 2, sd)/sqrt(n)
                 )

png("NumberOfDistinctEnhancers_vs_SampleSize.png", width = 6, height = 6, units = "in", res = 1000)
plotCI(x = dat$x, y = dat$y, li = dat$low, ui = dat$up, xlab = "Sample size", ylab = "Number of distinct enhancers")
dev.off()
