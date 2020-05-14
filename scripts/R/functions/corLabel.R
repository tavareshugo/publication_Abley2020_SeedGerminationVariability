#' Helper function to add correlation to graph
corLabel <- function(x, y, ci = FALSE, ...){
  cor_test <- cor.test(x, y, ...)

  # Get correlation, confidence interval and p-value
  r <- cor_test$estimate %>% prettyNum(nsmall = 2, digits = 2)
  r_ci <- cor_test$conf.int %>% prettyNum(nsmall = 2, digits = 2) %>% paste(collapse = ", ")
  r_p <- cor_test$p.value %>% prettyNum(digits = 2)

  if(ci){
    cor_label <- paste0("r = ", r, " [", r_ci, "]")
  } else {
    cor_label <- paste0("r = ", r)
  }
  return(cor_label)
}
