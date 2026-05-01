dmixsespc <- function(y, probs, mu, theta, logden = FALSE) {
  g <- length(probs)
  den <- matrix(nrow = dim(y)[1], ncol = g)
  for (j in 1:g)  den[, j] <- Directional::dsespc(y, mu[j, ], theta[j, ])
  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
