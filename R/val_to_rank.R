val_to_rank <- function(vals) {
  n <- length(vals)
  ranks <- 1:n
  return(ranks[order(vals)])
}
