np_order_stat <- function(size, alpha, delta) {
  violation_rates <- pbinom(q=0:(size-1), size=size, prob=1-alpha, lower.tail=F)
  return( which(violation_rates <= delta)[1] )
}
