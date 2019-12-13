library(recount)
library(MBCluster.Seq)
library(parallel)

set.seed(1)

options(mc.cores = snakemake@threads)

#download_study("SRP058719")
load("SRP058719/rse_gene.Rdata")

pheno <- colData(rse_gene)$title
names(pheno) <- colData(rse_gene)$run
pheno <- sapply(strsplit(pheno, " "), "[", 1)
pheno <- factor(pheno)
counts <- assay(scale_counts(rse_gene))

mb <- RNASeq.Data(counts, NULL, pheno[colnames(counts)])
c0 <- mclapply(c(2,3,4,5,6,7,8,9,10,12,15,20,50), function(k) lapply(1:5, function(i) KmeansPlus.RNASeq(mb, k, "nbinom")))
cls <- mclapply(c0, function(cc) {
  tmp <- lapply(cc, function(ci) {
    Cluster.RNASeq(mb, "nbinom", ci$centers)
  })
  logl <- sapply(tmp, function(x) lglk.cluster.one(mb, "nbinom", x$cluster))
  tmp[[which.max(logl)]]
})

results <- list(mb, cls)
saveRDS(results, snakemake@output[[1]])