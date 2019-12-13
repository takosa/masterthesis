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
p$foldchange <- "random"
jsonlite::write_json(p, output[[2]])

library(parallel)
devtools::load_all("benchDEwithClustering")
data(arab, package = "TCC")

all_results <- mclapply(seq_len(p$nsim), function(s) {
  set.seed(s)
  #####################################
  ### Generation of simulation data ###
  #####################################
  ## different levels and distributions of DE, where
  ## the fold-changes for DEGs are randomly sampled from 
  ## "1.2 + a gamma distribution with shape = 2.0 and scale = 0.5"
  ## mean fold-change is 2.2 (= 1.2 + 2.0 * 0.5)
  ## For obtaining Additinal file 1 Sheet 5
  #source("http://www.iu.a.u-tokyo.ac.jp/~kadota/TCC/TCC.simulation.R")
  fc.params <- data.frame(
    floor = c(1.2, 1.2, 1.2),
    shape = c(2.0, 2.0, 2.0),
    scale = c(0.5, 0.5, 0.5))
  fcm <- TCC::makeFCMatrix(Ngene = p$ngene, PDEG = p$pdeg, 
                           replicates = p$replicates, DEG.assign = p$assign,
                           fc.params = fc.params)
  tcc <- suppressWarnings(TCC::simulateReadCounts(Ngene = p$ngene,
                                                  PDEG = p$pdeg,
                                                  DEG.assign = p$assign,
                                                  DEG.foldchange = NULL,
                                                  replicates = p$replicates,
                                                  fc.matrix = fcm))
  
  # estimate size factors with EEE normalization in TCC package
  sizefactors <- get_sizefactors(counts = tcc$count, group = tcc$group$group, method = "eee")
  
  result <- list()
  result$simdata <- tcc
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 2\n")
  result$eee_mbcluster2 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 2, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 3\n")
  result$eee_mbcluster3 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 3, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 4\n")
  result$eee_mbcluster4 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 4, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 5\n")
  result$eee_mbcluster5 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 5, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 2\n")
  result$eee_mbcluster2_nstart5 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 2, nstart = 5, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 3\n")
  result$eee_mbcluster3_nstart5 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 3, nstart = 5, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 4\n")
  result$eee_mbcluster4_nstart5 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 4, nstart = 5, sizefactors = sizefactors)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 5\n")
  result$eee_mbcluster5_nstart5 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 5, nstart = 5, sizefactors = sizefactors)
 
  
   ## MBCluster.Seq
  #cat("start MBCluster.Seq, k = 6\n")
  #result$mbcluster6 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 6)
  ## MBCluster.Seq
  #cat("start MBCluster.Seq, k = 10\n")
  #result$mbcluster10 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 10)
  ## MBCluster.Seq
  #cat("start MBCluster.Seq, k = 20\n")
  #result$mbcluster20 <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 20)
  ## MBCluster.Seq with nstarts
  #cat("start MBCluster.Seq with nstart\n")
  #result$mbcluster3_nstart <- run_mbcluster(counts = tcc$count, group = tcc$group$group, k = 4, nstart = 5)
  ## edgeR
  #cat("start edgeR\n")
  #result$edger <- run_edegr(counts = tcc$count, group = tcc$group$group, coef = 2)
  ## DESeq2
  #cat("start DESeq2\n")
  #result$deseq2 <- run_deseq2(counts = tcc$count, group = tcc$group$group, test = "LRT", reduced = ~1)
  ## TCC
  #cat("start TCC\n")
  #result$tcc <- run_tcc(counts = tcc$count, group = tcc$group$group, X = "tmm", Y = "edger", Z = "edger", iter = 3)
  ## DREAMSeq
  #cat("start DREAMSeq\n")
  #result$dreamseq <- run_dreamseq(counts = tcc$count, group = tcc$group$group)
  
  result
})

rm(arab)
saveRDS(all_results, file = output[[1]])
