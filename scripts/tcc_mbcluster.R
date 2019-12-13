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
  # estimate size factors with EEE normalization in TCC package
  sizefactors <- get_sizefactors(counts = tcc$count, group = tcc$group$group, method = "eee")
  
  # EEE-MBCluster.Seq
  cat("start MBCluster.Seq, k = 3, with EEE normalization\n")
  result$eee_mbcluster3 <- run_mbcluster(counts = tcc$count, group = tcc$group$group,
                                     sizefactors = sizefactors, k = 3)
  # EEE-MBCluster.Seq, nstart = 5
  cat("start MBCluster.Seq, k = 3, nstart = 5, with EEE normalization\n")
  result$eee_mbcluster3_nstart <- run_mbcluster(counts = tcc$count, group = tcc$group$group,
                                     sizefactors = sizefactors, k = 3, nstart = 5)
  
  cat("start MBCluster.Seq, k = 3, nstart = 5, with EEE normalization\n")
  result$zero_mbcluster3 <- run_mbcluster(counts = tcc$count, group = tcc$group$group,
                                          sizefactors = rep(1, ncol(tcc$count)), k = 3)
  
  result
})

rm(arab)
saveRDS(all_results, file = output[[1]])
