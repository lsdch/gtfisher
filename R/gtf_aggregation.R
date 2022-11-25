#' @import dplyr

#' @export
fisher_score = function(pvals) {
  -2 * sum(log(pvals))
}

gtf_pval = function(gtf_score, k, size, n, pb = function() {
                      NULL
                    }) {
  # Workaround for group having less than `k` p-values
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


#' Compute aggregated p-value per group of p-values resulting from independent tests.
#'
#' Computation of aggregated p-values is performed using the so-called
#' "Groupwise Truncated Fisher" method.
#' It consists in compute a Fisher score on the best ranking p-values within each group,
#' then use an empirical distribution of this score to calculate a p-value
#' for the corresponding group.
#' @param dataset a tibble or dataframe
#' @param group_col column of the grouping variable
#' @param pval_col column of p-values
#' @param ... filtering condition as in dplyr::filter
#' @param k the number of best ranking p-values to consider in each group
#' @param n the sample size to use in the empirical null distribution
#' @examples
#' # no progress bar
#' group_pvals = gtf_predict(pvals, group, p)
#'
#' # with progress bar
#' handlers("progress") # pick your favorite one, see `progressr` vignette
#' with_progress({
#'   group_pvals = gtf_predict(pvals, group, p)
#' })
#'
#' @export
gtf_predict = function(dataset, group_col, pval_col, ..., k = 5, n = 10000) {
  grp_data = dataset %>% group_by({{ group_col }})

  if (requireNamespace("progress", quietly = TRUE)) {
    p <- progressr::progressor(steps = n_groups(grp_data))
  } else {
    p = function() NULL
  }

  grp_data %>%
    filter(...) %>%
    summarize(
      size = n(),
      gtf_score = fisher_score(DescTools::Small({{ pval_col }}, k)),
      gtf_pval = gtf_pval(gtf_score, k = k, size = size, n = n, pb = p),
      .groups = "keep"
    ) %>%
    arrange(gtf_pval, desc(gtf_score))
}
