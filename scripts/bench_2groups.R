#snakemake@input           #: list()
output <- snakemake@output          #: list()
p <- snakemake@params[[1]]          #: list()
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
#p <- list(nsim = 2,
#          ngene = 10000,
#          pdeg = 0.05,
#          assign = c(1/2, 1/2),
#          foldchange = c(4, 4),
#          replicates = c(3, 3),
#          group = NULL,
#          fcmatrix = NULL)
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
jsonlite::write_json(p, output[[2]])

library(parallel)
devtools::load_all("benchDEwithClustering")
data(arab, package = "TCC")

all_results <- mclapply(seq_len(p$nsim), function(s) {
  set.seed(s)
  tcc <- suppressWarnings(TCC::simulateReadCounts(Ngene = p$ngene,
                                                  PDEG = p$pdeg,
                                                  DEG.assign = p$assign,
                                                  DEG.foldchange = p$foldchange,
                                                  replicates = p$replicates,
                                                  group = p$group,
                                                  fc.matrix = p$fcmatrix))
  
  result <- list()
  result$simdata <- tcc
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 2\n")
  result$mbcluster2 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 3\n")
  result$mbcluster3 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 3)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 4\n")
  result$mbcluster4 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 4)
  # MBCluster.Seq with nstarts
  cat("start MBCluster.Seq with nstart\n")
  result$mbcluster_nstart <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 3, nstart = 5)
  # edgeR
  cat("start edgeR\n")
  result$edger <- run_edger(counts = tcc$count, group = tcc$group$group, coef = 2)
  # DESeq2
  cat("start DESeq2\n")
  result$deseq2 <- run_deseq2(counts = tcc$count, group = tcc$group$group, test = "LRT", reduced = ~1)
  # TCC
  cat("start TCC\n")
  result$tcc <- run_tcc(counts = tcc$count, group = tcc$group$group, X = "tmm", Y = "edger", Z = "edger", iter = 3)
  # DREAMSeq
  cat("start DREAMSeq\n")
  result$dreamseq <- run_dreamseq(counts = tcc$count, group = tcc$group$group)
  
  result
})

rm(arab)
saveRDS(all_results, file = output[[1]])
