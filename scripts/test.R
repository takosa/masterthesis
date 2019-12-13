library(TCC)
library(MBCluster.Seq)

tcc <- simulateReadCounts(20000, 0.2, c(.5, .5), c(4, 4), c(3, 3))
mb <- RNASeq.Data(tcc$count, rep(0, 6), gl(2,3))
c0 <- lapply(1:10, function(x) KmeansPlus.RNASeq(mb, 10, "nbinom"))
cls <- lapply(c0, function(x) Cluster.RNASeq(mb, "nbinom", x$centers))
logL <- sapply(cls, function(x) lglk.cluster(mb, "nbinom", x$cluster))
plot(logL, type = "b", xlab = "iteration", ylab = "log likelihood")
for (i in 1:10) {
  plot(cls[[i]]$centers, asp = 1, col = 2)
  abline(0, -1, lty = 2)
  abline(0, 0, lty = 2)
  abline(v = 0, lty = 2)
  points(z[1], z[2], pch = 4, cex = 3, col = 3)
  points(z[2], z[1], pch = 4, cex = 3, col = 3)
  points(0, 0, pch = 4, cex = 3, col = 3)
  readline()
}

x <- colsum(tcc$count, c(1,1,1,2,2,2))
x <- cbind(log(x[, 1]) + log(x[, 2]), log(x[, 1]) - log(x[, 2]))
col <- cls[[which.max(logL)]]$cluster
mycol <- c("#ff634799", "#4169e199", "#3cb37199", "#ffd70099", "#2f4f4f99",
           "#ff00ff99", "#9400d399", "#00ffff99", "#6b8e2399", "#4b008299",
           "#fa807255")[col]
plot(x, col = mycol, pch = 16, cex = 0.7)

tcc <- simulateReadCounts(20000, 0.2, c(.5, .5), c(4, 4), c(3, 3))
mb <- RNASeq.Data(tcc$count, rep(0, 6), gl(2,3))
c0 <- lapply(1:10, function(x) KmeansPlus.RNASeq(mb, 3, "nbinom"))
cls <- lapply(c0, function(x) Cluster.RNASeq(mb, "nbinom", x$centers))
logL <- sapply(cls, function(x) lglk.cluster(mb, "nbinom", x$cluster))
plot(logL, type = "b", xlab = "iteration", ylab = "log likelihood")
for (i in 1:10) {
  plot(cls[[i]]$centers, asp = 1, col = 2)
  abline(0, -1, lty = 2)
  abline(0, 0, lty = 2)
  abline(v = 0, lty = 2)
  points(z[1], z[2], pch = 4, cex = 3, col = 3)
  points(z[2], z[1], pch = 4, cex = 3, col = 3)
  points(0, 0, pch = 4, cex = 3, col = 3)
  readline()
}

x <- colsum(tcc$count, c(1,1,1,2,2,2))
x <- cbind(log(x[, 1]) + log(x[, 2]), log(x[, 1]) - log(x[, 2]))
col <- cls[[which.max(logL)]]$cluster
mycol <- c("#ff634799", "#4169e199", "#3cb37199", "#ffd70099", "#2f4f4f99",
           "#ff00ff99", "#9400d399", "#00ffff99", "#6b8e2399", "#4b008299",
           "#fa807255")[col]
plot(x, col = mycol, pch = 16, cex = 0.7)

load("~/Downloads/rse_gene.Rdata")
rse_gene
colData(rse_gene)$geo_accession

dplyr::select(tibble::as_tibble(colData(rse_gene)), sample, experiment, run,
              mapped_read_count, auc, title, characteristics)

library(recount)
geo_characteristics(colData(rse_gene))
table(geo_characteristics(colData(rse_gene)))


