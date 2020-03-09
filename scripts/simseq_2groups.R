#snakemake@input           #: list()
output <- snakemake@output          #: list()
#p <- snakemake@params[[1]]          #: list()
#snakemake@wildcards       #: list()
threads <- snakemake@threads         #: num 1
log <- snakemake@log             #: list()
#snakemake@resources       #: list()
#snakemake@config          #: list()
#snakemake@rule            #: chr "all"
#snakemake@bench_iteration #: num NA
#snakemake@scriptdir       #: chr "/scripts"
#snakemake@source          #:function (...)

## if run this script directly without snakemake
#output <- list("out/out.rds")
#threads <- 1
#log <- list()

if (length(log) > 0) {
  if (!dir.exists(dname <- dirname(log[[1]])))
    dir.create(dname, recursive = TURE)
  if (length(log) == 2) {
    zz1 <- file(log[[1]], open = "wt")
    sink(zz1)
    zz1 <- file(log[[2]], open = "wt")
    sink(zz1, type = "message")
  } else{
    zz <- file(log[[1]], open = "wt")
    sink(zz)
    sink(zz, type = "message")
  }
}

options(mc.cores = threads)

library(parallel)
devtools::load_all("benchDEwithClustering")
data(kidney, package = "SimSeq")
nf <- apply(kidney$counts, 2, quantile, 0.75)
sorted <- SimSeq::SortData(counts = kidney$counts, replic = kidney$replic,
                              treatment = kidney$treatment, sort.method = "paired",
                              norm.factors = nf)
probs <- SimSeq::CalcPvalWilcox(counts = sorted$counts, treatment = sorted$treatment,
                                replic = sorted$replic, sort.method = "paired",
                                sorted = TRUE, norm.factors = nf, exact = FALSE)
weights <- 1 - fdrtool::fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr

all_results <- mclapply(seq_len(12), function(s) {
  set.seed(s)
  sim <- SimSeq::SimData(counts = sorted$counts, replic = sorted$replic,
                         treatment = sorted$treatment, sort.method = "paired",
                         k.ind = 5, n.genes = 10000,n.diff = 1000, weights = weights,
                         norm.factors = nf)
  
  result <- list()
  result$simdata <- sim
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 7\n")
  result$mbcluster7 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 7)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 8\n")
  result$mbcluster8 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 8)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 9\n")
  result$mbcluster9 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 9)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 10\n")
  result$mbcluster10 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 10)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 11\n")
  result$mbcluster11 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 11)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 12\n")
  result$mbcluster12 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 12)
  
  # edgeR dispersion
  disp <- get_dispersions(counts = sim$counts, group = as.factor(sim$treatment), method = "edger")
  size <- get_sizefactors(counts = sim$counts, group = as.factor(sim$treatment), method = "tmm")
  
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 7, edgeR disp\n")
  result$mbcluster_disp7 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 7, nstart = 5, dispersions = disp, sizefactors = size)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 8, edgeR disp\n")
  result$mbcluster_disp8 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 8, nstart = 5, dispersions = disp, sizefactors = size)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 9, edgeR disp\n")
  result$mbcluster_disp9 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 9, nstart = 5, dispersions = disp, sizefactors = size)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 10, edgeR disp\n")
  result$mbcluster_disp10 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 10, nstart = 5, dispersions = disp, sizefactors = size)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 11, edgeR disp\n")
  result$mbcluster_disp11 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 11, nstart = 5, dispersions = disp, sizefactors = size)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 12, edgeR disp\n")
  result$mbcluster_disp12 <- run_mbcluster(counts = sim$counts, group = as.factor(sim$treatment), k = 12, nstart = 5, dispersions = disp, sizefactors = size)
  
  # edgeR
  cat("start edgeR\n")
  result$edger <- run_edger(counts = sim$counts, group = as.factor(sim$treatment), coef = 2)
  # DESeq2
  cat("start DESeq2\n")
  result$deseq2 <- run_deseq2(counts = sim$counts, group = as.factor(sim$treatment), test = "LRT", reduced = ~1)
  # TCC
  cat("start TCC\n")
  result$tcc <- run_tcc(counts = sim$counts, group = as.factor(sim$treatment), X = "tmm", Y = "edger", Z = "edger", iter = 3)
  # DREAMSeq
  cat("start DREAMSeq\n")
  result$dreamseq <- run_dreamseq(counts = sim$counts, group = as.factor(sim$treatment))
  
  result
})

saveRDS(all_results, file = output[[1]])
