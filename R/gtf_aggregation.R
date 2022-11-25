#' @import dplyr

#' @export
fisher_score = function(pvals) {
  -2 * sum(log(pvals))
}

gtf_pval = function(gtf_score, k, size, n, pb = function() {
                      NULL
                    }) {
  # Workaround for genes having less than `k` non constant sites
  k = min(k, size)
  # Simulated null distribution
  null_scores = sapply(c(1:n), function(x) {
    # draw the first k minimum values from uniform
    # https://en.wikipedia.org/wiki/Order_statistic#Probability_distributions_of_order_statistics
    rbeta(k, c(1:k), size + 1 - c(1:k)) %>% fisher_score()
  })
  # Report progress
  pb()
  # Test for deviation from the null
  plnorm(gtf_score,
    meanlog = mean(log(null_scores)),
    sdlog = sd(log(null_scores)),
    lower.tail = F
  )
}


#' Compute gene-level p-values from Pelican site-wise p-values.
#'
#' Computation of gene p-values is performed using the so-called
#' "Genewise Truncated Fisher" method.
#' It consists in compute a Fisher score on the best ranking p-values within each gene,
#' and use an empirical distribution of this score to calculate a p-value
#' for the corresponding gene.
#' @param dataset a tibble or dataframe that contains site p-values for a collection of genes
#' @param ali_col column of gene identifier
#' @param pval_col column of site p-values
#' @param naa_col column that contains the number of distinct amino acid at each site
#' @param k the number of best ranking p-values to consider at each gene
#' @param n the sample size to use in the empirical null distribution
#' @examples
#' site_pvals = read_tsv("/path/to/pelican/output/all_sites.tsv")
#' # no progress bar
#' gene_predictions = site_pvals %>% gtf_predict()
#'
#' # with progress bar
#' handlers("progress") # pick your favorite one, see `progressr` vignette
#' with_progress({
#'   gene_predictions = site_pvals %>% gtf_predict()
#' })
#'
#' @export
gtf_predict = function(dataset,
                       ali_col = alignment,
                       pval_col = aagtr_pval,
                       naa_col = naa,
                       k = 5,
                       n = 10000) {
  grp_data = dataset %>% group_by({{ ali_col }})

  if (requireNamespace("progress", quietly = TRUE)) {
    p <- progressr::progressor(steps = n_groups(grp_data))
  } else {
    p = function() NULL
  }

  grp_data %>%
    summarize(
      size_non_constant = sum({{ naa_col }} > 1),
      gtf_score = fisher_score(DescTools::Small(aagtr_pval, k)),
      gtf_pval = gtf_pval(gtf_score, k, size_non_constant, n = n, pb = p),
      .groups = "keep"
    ) %>%
    arrange(gtf_pval, desc(gtf_score))
}
