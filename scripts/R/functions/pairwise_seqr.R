#' Run QTLseqr analysis for all pairwise comparisons of samples
#'
#' @param file the file with counts
#' @param windowSize the window size parameter passed to QTLseqr functions
#' @param replications the number of replicates for the threshold estimation
#' @param pool_sizes a named vector with size of each pool. For example `c(pool1 = 100, pool2 = 20)`
#' @param ... other parameters passed to QTLseqr functions
pairwise_seqr <- function(file,
                          windowSize = 1e6,
                          replications = 1e4,
                          filter = 0,
                          pool_sizes, ...){

  require(dplyr); require(readr); require(purrr); require(stringr)

  # read counts
  counts <- read_tsv(file)

  # Extract sample names from table columns
  samples <- counts %>%
    select(matches("\\.AD")) %>%
    colnames() %>%
    str_remove("\\.AD")

  message("The following samples were detected: ", paste(samples, collapse = ", "))

  # Make pairwise combinations
  comparisons <- combn(samples, 2) %>% as.data.frame %>% as.list %>% map(as.character)

  # Run QTLseqr pipeline for each pairwise comparison
  bsa_qtl <- map_df(comparisons, function(i){

    cat("Pools: ", pool_sizes[i], "\n")
    cat(i, "\ns")

    # Run QTLseqr pipeline
    importFromGATK(file, lowBulk = i[1], highBulk = i[2]) %>%
      runGprimeAnalysis(windowSize = windowSize) %>%
      runQTLseqAnalysis(windowSize = windowSize, replications = replications,
                        bulkSize = pool_sizes[i], filter = filter, ...) %>%
      mutate(comparison = paste0(i[2], " - ", i[1]))
  })

  return(bsa_qtl)
}
