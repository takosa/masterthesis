input <- snakemake@input           #: list()
output_file <- snakemake@output[[1]]          #: list()
#p <- snakemake@params          #: list()
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

result_input <- input[[1]]
param_input <- input[[2]]

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
param <- jsonlite::read_json(param_input)
all_results <- readRDS(result_input)

library(dplyr)
library(purrr)
library(yardstick)

result <- all_results %>%
  imap(function(res, i) {
    simdata <- data.frame(truth = factor(res$simdata$simulation$trueDEG, levels = c(1, 0))) %>%
      tibble::rownames_to_column("gene_id") %>% as_tibble()
    imap(res[-1], function(re, method) {
      re$summary %>%
        right_join(simdata, by = "gene_id") %>%
        mutate(Method = method)
    }) %>% bind_rows() -> tmp
    mutate(tmp, Iter = as.integer(i))
  }) %>% bind_rows()

Ngene <- flatten_int(param$ngene)
PDEG <- flatten_dbl(param$pdeg)
Assign <- stringr::str_flatten(map_chr(param$assign, ~sprintf("%.2f", .)), collapse = "_")
Foldchange <- as.character(param$foldchange[[1]])
Replicates <- param$replicates[[1]]
result <- mutate(result,
                 Ngene = Ngene,
                 PDEG = PDEG,
                 Assign = Assign,
                 Foldchange = Foldchange,
                 Replicates = Replicates) %>%
  as_tibble()

multi_metric <- metric_set(sens, spec, recall, precision, mcc, j_index, 
                           f_meas, accuracy, kap, ppv, npv, bal_accuracy,
                           detection_prevalence, roc_auc, pr_auc, 
                           average_precision, gain_capture, mn_log_loss)

met <- result %>%
  mutate(prediction = factor(if_else(q.value <= 0.05, 1, 0), levels = c(1, 0))) %>%
  mutate(p.value.inv = 1.0 - p.value) %>%
  group_by(Ngene, PDEG, Assign, Foldchange, Replicates, Method, Iter) %>%
  multi_metric(truth = truth, estimate = prediction, p.value.inv)

tmp <- list(result = result, metrics = met)
saveRDS(tmp, output_file)


