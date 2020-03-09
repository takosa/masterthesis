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
  genes.de <- sample(seq_len(nrow(kidney$counts)), size = 1000, prob = weights) # Sample all DE genes
  DE1 <- genes.de[1:333] # Sample DE genes with first trt diff
  DE2 <- genes.de[334:666] # Sample DE genes with sec trt diff
  DE3 <- genes.de[667:1000] # Sample DE genes with third trt diff
  EE <- sample( seq_len(nrow(kidney$counts))[-genes.de], size = 4000) #Sample EE genes
  genes.tot <- c(EE, genes.de)
  genes.de1 <- union(DE2, EE) #Assign DE genes for first sim
  genes.de2 <- union(DE2, DE3) #Assign DE genes for second sim
  data.sim1 <- SimSeq::SimData(counts = kidney$counts, replic = kidney$replic,
                               treatment = kidney$treatment, sort.method = "paired",
                               k.ind = 5L, genes.select = genes.tot,
                               genes.diff = genes.de1, weights = weights, norm.factors = nf)
  #remove pairs of columns used in first simulation
  cols.rm <- c(data.sim1$col[1:(2*5L)], data.sim1$col[1:(2*5L)] + 1)
  counts.new <- kidney$counts[, -cols.rm]
  nf.new <- nf[-cols.rm]
  replic.new <- kidney$replic[-cols.rm]
  treatment.new <- kidney$treatment[-cols.rm]
  ### Set switch.trt = TRUE for second sim
  data.sim2 <- SimSeq::SimData(counts = counts.new, replic = replic.new, treatment = treatment.new,
                               sort.method = "paired", k.ind = 5L, genes.select = genes.tot,
                               genes.diff = genes.de2, weights = weights, norm.factors = nf.new,
                               switch.trt = TRUE)
  ### Remove first k.ind entries from first sim and combine two count matrices
  counts.sim <- cbind(data.sim1$counts[, -(1:5L)], data.sim2$counts)
  ### treatment group levels for simulated matrix
  trt.grp <- rep(NA, 5000)
  trt.grp[is.element(data.sim1$genes.subset, DE1)] <- "DE_First_Trt"
  trt.grp[is.element(data.sim1$genes.subset, DE2)] <- "DE_Second_Trt"
  trt.grp[is.element(data.sim1$genes.subset, DE3)] <- "DE_Third_Trt"
  trt.grp[is.element(data.sim1$genes.subset, EE)] <- "EE"
  treat <- gl(3, 5)
  result <- list()
  result$simdata <- list(counts = counts.sim, truth = trt.grp)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 5\n")
  result$mbcluster5 <- run_mbcluster(counts = counts.sim, group = treat, k = 5)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 7\n")
  result$mbcluster7 <- run_mbcluster(counts = counts.sim, group = treat, k = 7)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 11\n")
  result$mbcluster11 <- run_mbcluster(counts = counts.sim, group = treat, k = 11)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 15\n")
  result$mbcluster15 <- run_mbcluster(counts = counts.sim, group = treat, k = 15)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 21\n")
  result$mbcluster21 <- run_mbcluster(counts = counts.sim, group = treat, k = 21)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 29\n")
  result$mbcluster29 <- run_mbcluster(counts = counts.sim, group = treat, k = 29)
  
  # edgeR dispersion
  disp <- get_dispersions(counts = counts.sim, group = treat, method = "deseq2")
  size <- get_sizefactors(counts = counts.sim, group = treat, method = "tmm")
  
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 5, edgeR disp\n")
  result$mbcluster_disp5 <- run_mbcluster(counts = counts.sim, group = treat, k = 5, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 7, edgeR disp\n")
  result$mbcluster_disp7 <- run_mbcluster(counts = counts.sim, group = treat, k = 7, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 11, edgeR disp\n")
  result$mbcluster_disp11 <- run_mbcluster(counts = counts.sim, group = treat, k = 11, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 15, edgeR disp\n")
  result$mbcluster_disp15 <- run_mbcluster(counts = counts.sim, group = treat, k = 15, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 21, edgeR disp\n")
  result$mbcluster_disp21 <- run_mbcluster(counts = counts.sim, group = treat, k = 21, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  # MBCluster.Seq
  cat("start MBCluster.Seq, k = 29, edgeR disp\n")
  result$mbcluster_disp29 <- run_mbcluster(counts = counts.sim, group = treat, k = 29, nstart = 5, dispersions = disp, sizefactors = size, select_type = 2)
  
  # edgeR
  cat("start edgeR\n")
  result$edger <- run_edger(counts = counts.sim, group = treat, coef = 2:3)
  # DESeq2
  cat("start DESeq2\n")
  result$deseq2 <- run_deseq2(counts = counts.sim, group = treat, test = "LRT", reduced = ~1)
  
  result
})

saveRDS(all_results, file = output[[1]])
