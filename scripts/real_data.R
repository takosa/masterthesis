library(recount)
library(MBCluster.Seq)
library(parallel)

devtools::load_all("benchDEwithClustering")
set.seed(1)

options(mc.cores = snakemake@threads)

#download_study("SRP058719")
load("SRP058719/rse_gene.Rdata")

pheno <- colData(rse_gene)$title
names(pheno) <- colData(rse_gene)$run
pheno <- sapply(strsplit(pheno, " "), "[", 1)
pheno <- factor(pheno)
counts <- assay(scale_counts(rse_gene))

K <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 50)
res <- mclapply(K, function(k) run_mbcluster(counts, pheno, k = k, nstart = 5), mc.set.seed = 1L)
names(res) <- paste0("mbcluster", K)

res$deseq2 <- run_deseq2(counts, pheno, test = "LRT", reduced = ~1)
res$edger <- run_edger(counts, pheno, coef = 2:10)


saveRDS(res, snakemake@output[[1]])

