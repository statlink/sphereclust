dmixkent <- function(y, probs, G, param, logden = FALSE) {
  g <- length(probs)
  den <- matrix(nrow = dim(y)[1], ncol = g)
  for (j in 1:g)  den[, j] <- Directional::dkent(y, G[[ j ]], param[j, ])
  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
