dmixesag <- function(y, probs, mu, gam, logden = FALSE) {
  g <- length(probs)
  den <- matrix(nrow = dim(y)[1], ncol = g)
  for (j in 1:g)  den[, j] <- Directional::desag(y, mu[j, ], gam[j, ])
  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
