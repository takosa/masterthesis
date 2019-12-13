#' @import magrittr
run_mbcluster <- function(counts, group, sizefactors = NULL, k = 3, nstart = 1, dispersions = NULL, select_type = 1) {

  ptm <- proc.time()
  ###
  if (is.null(sizefactors)) {
    mb <- MBCluster.Seq::RNASeq.Data(Count = counts, Normalizer = NULL, Treatment = group)
  } else {
    mb <- MBCluster.Seq::RNASeq.Data(Count = counts, Normalizer = log(sizefactors), Treatment = group)
  }
  c0 <- lapply(seq_len(nstart), function(i) {
    MBCluster.Seq::KmeansPlus.RNASeq(data = mb, nK = k, model = "nbinom")
  })
  cls <- lapply(c0, function(cc) {
    MBCluster.Seq::Cluster.RNASeq(data = mb, model = "nbinom", centers = cc$centers,
                                  method = "EM", iter.max = 50, TMP = NULL)
  })
  logL <- sapply(cls, function(cc) {
    MBCluster.Seq::lglk.cluster(data = mb, model = "nbinom", cluster = cc$cluster)
  })
  ###
  ptm <- proc.time() - ptm
  if (select_type == 1) {
    res <- cls[[which.max(logL)]]
  } else if (select_type == 2) {
    res <- cls[[which.min(sapply(cls, function(x) min(rowSums(x$centers^2))))]]
  }
  nde_idx <- which.min(rowSums(res$centers ^ 2))
  pval <- res$probability[, nde_idx]
  summary_res <- data.frame(gene_id = rownames(counts),
                            p.value = pval,
                            stringsAsFactors = FALSE) %>%
    dplyr::arrange(p.value) %>%
    dplyr::mutate(q.value = dplyr::cummean(p.value))

  rtn <- list()
  rtn$time <- ptm
  rtn$result <- list(data = mb, cluster = res, log_likelihood = max(logL))
  rtn$summary <- summary_res

  rtn
}

#' @import magrittr
run_edger <- function(counts, group, coef) {

  ptm <- proc.time()
  ###
  y <- edgeR::DGEList(counts = counts, group = group)
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~group)
  y <- edgeR::estimateDisp(y, design)
  # use glmQLFit rather than glmFit because glmQLFit is recommended in edgeR user guid
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = coef)
  res <- edgeR::topTags(qlf, n = nrow(counts), sort.by = "none")
  #fit <- glmFit(y, design)
  #lrt <- glmLRT(fit, coef = coef)
  #edgeR::topTags(lrt)
  ###
  ptm <- proc.time() - ptm

  summary_res <- res$table %>% tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, p.value = PValue, q.value = FDR)

  rtn <- list()
  rtn$time <- ptm
  rtn$result <- res
  rtn$summary <- summary_res

  rtn
}

#' @import magrittr
#' @importClassesFrom TCC TCC
run_tcc <- function(counts, group, coef, X, Y, Z, iter) {

  ptm <- proc.time()
  ###
  tcc <- new("TCC", count = counts, group = group)
  tcc <- TCC::calcNormFactors(tcc, norm.method = X, test.method = Y, iter = iter)
  tcc <- TCC::estimateDE(tcc, test.method = Z)
  res <- TCC::getResult(tcc)
  ###
  ptm <- proc.time() - ptm

  summary_res <- res %>% dplyr::select(gene_id, p.value, q.value) %>%
    dplyr::mutate(gene_id = as.character(gene_id))

  rtn <- list()
  rtn$time <- ptm
  rtn$result <- res
  rtn$summary <- summary_res

  rtn
}

#' @import magrittr
run_deseq2 <- function(counts, group, test, reduced) {

  ptm <- proc.time()
  ###
  coldata <- data.frame(group = group)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~group)
  dds <- DESeq2::DESeq(dds, test = test, reduced = reduced)
  res <- DESeq2::results(dds)
  ###
  ptm <- proc.time() - ptm

  summary_res <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, p.value = pvalue, q.value = padj)
  rtn <- list()
  rtn$time <- ptm
  rtn$result <- res
  rtn$summary <- summary_res

  rtn
}

#' @import magrittr
run_dreamseq <- function(counts, group) {

  ptm <- proc.time()
  ###
  cnt <- DREAMSeq::DEGCountSet(counts, group)
  cnt <- DREAMSeq::DREAMSeq(cnt, model = "DP", parallel = F)
  res <- DREAMSeq::exResult(cnt)
  ###
  ptm <- proc.time() - ptm

  summary_res <- res %>% dplyr::select(gene_id = geneId, p.value = pValue, q.value = pAdj)
  rtn <- list()
  rtn$time <- ptm
  rtn$result <- res
  rtn$summary <- summary_res

  rtn
}

#' @import magrittr
run_bayseq <- function(counts, group, model = NULL) {

  ptm <- proc.time()
  ###
  cd <- new("countData", data = counts, replicates = group)
  if (is.null(model)) cd <- baySeq::allModels(cd)
  else baySeq::groups(cd) <- model
  baySeq::libsizes(cd) <- baySeq::getLibsizes(cd, estimationType = "edgeR")
  cd@annotation <- data.frame(gene_id = rownames(counts))
  cd <- baySeq::getPriors.NB(cd, samplesize = 3000, estimation = "QL", cl = NULL)
  cd <- baySeq::getLikelihoods(cd, cl = NULL, bootStraps = 3, verbose = FALSE)
  ###
  ptm <- proc.time() - ptm

  res <- cd
  summary_res <- res <- topCounts(cd, group = 1)

  summary_res <- res %>% dplyr::select(gene_id = geneId, p.value = pValue, q.value = pAdj)
  rtn <- list()
  rtn$time <- ptm
  rtn$result <- res

  rtn
}

mrnFactors <- function(rawCounts,conditions) {
  rawCounts <- as.matrix(rawCounts)
  totalCounts <- colSums(rawCounts)
  normFactors <- totalCounts
  medianRatios <- rep(1,length(conditions))
  names(medianRatios) <- names(normFactors)
  if (sum(conditions==1)>1)
    meanA <- apply(rawCounts[,conditions==1]%*%diag(1/totalCounts[conditions==1]),1,mean)
  else
    meanA <- rawCounts[,conditions==1]/totalCounts[conditions==1]
  for (i in 2:max(conditions)) {
    if (sum(conditions==i)>1)
      meanB <- apply(rawCounts[,conditions==i]%*%diag(1/totalCounts[conditions==i]),1,mean)
    else
      meanB <- rawCounts[,conditions==i]/totalCounts[conditions==i]
    meanANot0 <- meanA[meanA>0&meanB>0]
    meanBNot0 <- meanB[meanA>0&meanB>0]
    ratios <- meanBNot0/meanANot0
    medianRatios[conditions==i] <- median(ratios)
    normFactors[conditions==i] <- medianRatios[conditions==i]*totalCounts[conditions==i]
  }
  medianRatios <- medianRatios/exp(mean(log(medianRatios)))
  normFactors <- normFactors/exp(mean(log(normFactors)))
  return(list(medianRatios=medianRatios,normFactors=normFactors))
}

get_sizefactors <- function(counts, group, method = c("tmm", "rle", "mrn", "eee")) {

  method <- match.arg(method)
  switch (method,
    tmm = {
      normfactors <- edgeR::calcNormFactors(counts)
      libsizes <- normfactors * colSums(counts)
      sizefactors <- libsizes / mean(libsizes)
    },
    rle = {
      sizefactors <- DESeq2::estimateSizeFactorsForMatrix(counts = counts)
    },
    mrn = {
      group <- as.integer(group)
      sizefactors <- mrnFactors(rawCounts = counts, conditions = group)$normFactors
    },
    eee = {
      tcc <- new("TCC", count = counts, group = group)
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iter = 3)
      libsizes <- tcc$norm.factors * colSums(counts)
      sizefactors <- libsizes / mean(libsizes)
    }
  )
  sizefactors
}


