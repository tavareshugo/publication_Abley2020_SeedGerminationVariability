#' Plot linear model diagnostics for a marker in the QTL scan
#'
#' @param x an object of class MagicGenPhen (from MagicHelpR)
#' @param phen the phenotype to fit
#' @param marker the marker to fit the model to
qtlDiag <- function(x, phen, marker){
  if(length(marker) != 1 | length(phen) != 1) stop("One phenotype and marker need to be given.")
  
  P <- getPhenotypes(x)[[phen]]
  G <- getGenotypes(x)[[marker]]
  
  par(mfrow = c(2, 2))
  lm(P ~ G) %>% plot()
}
