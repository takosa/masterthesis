threads <- snakemake@threads
options(mc.cores = threads)

log <- snakemake@log
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


library(tidyverse)

summarize_result <- read_rds(snakemake@input[[1]])$result

res <- list()

res$q <- summarize_result %>%
  select(-Ngene, -PDEG, -Assign, -Foldchange, -Replicates, -p.value) %>%
  nest(data = c(gene_id, q.value, truth)) %>% 
  mutate(data = map(data, function(d) {
    x <- seq(0, 1, length.out = 100)
    d <- arrange(d, q.value) %>%
      mutate(fdr = cummean(as.integer(truth == 0)))
    y <- approxfun(d$q.value, d$fdr, rule = 2, ties = "mean")(x)
    tibble(x = x, y = y)
  })) %>%
  unnest() %>%
  nest(data = c(Iter, x, y)) %>%
  mutate(data = map(data, function(d) {
    d %>% group_by(x) %>% summarise(y = mean(y))
  })) %>%
  unnest() %>%
  mutate(Ngene = as.character(unique(summarize_result$Ngene)),
         PDEG  = as.character(unique(summarize_result$PDEG)),
         Assign  = as.character(unique(summarize_result$Assign)),
         Foldchange  = as.character(unique(summarize_result$Foldchange)),
         Replicates  = as.character(unique(summarize_result$Replicates)))


res$p <- summarize_result %>%
  select(-Ngene, -PDEG, -Assign, -Foldchange, -Replicates, -q.value) %>%
  nest(data = c(gene_id, p.value, truth)) %>% 
  mutate(data = map(data, function(d) {
    x <- seq(0, 1, length.out = 100)
    d <- arrange(d, p.value) %>%
      mutate(fdr = cummean(as.integer(truth == 0)))
    y <- approxfun(d$p.value, d$fdr, rule = 2, ties = "mean")(x)
    tibble(x = x, y = y)
  })) %>%
  unnest() %>%
  nest(data = c(Iter, x, y)) %>%
  mutate(data = map(data, function(d) {
    d %>% group_by(x) %>% summarise(y = mean(y))
  })) %>%
  unnest() %>%
  mutate(Ngene = as.character(unique(summarize_result$Ngene)),
         PDEG  = as.character(unique(summarize_result$PDEG)),
         Assign  = as.character(unique(summarize_result$Assign)),
         Foldchange  = as.character(unique(summarize_result$Foldchange)),
         Replicates  = as.character(unique(summarize_result$Replicates)))


write_rds(res, snakemake@output[[1]])
