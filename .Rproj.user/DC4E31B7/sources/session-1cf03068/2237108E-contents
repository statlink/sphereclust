rmixsespc <- function(n, probs, mu, theta) {
  ## n is the sample size
  ## prob is a vector with the mixing probabilities
  ## mu is a matrix with with the mean directions
  ## k is a vector with the concentration parameters
  p2 <- c( 0, cumsum(probs) )
  p <- ncol(mu)  ## dimensionality of the data
  u <- rangen::Runif(n)
  g <- length(probs)  ## how many clusters are there
  ina <- as.numeric( cut(u, breaks = p2) )  ## the cluster of each observation
  ina <- sort(ina)
  nu <- tabulate(ina)  ## frequency table of each cluster
  y <- array( dim = c(n, p, g) )
  for (j in 1:g)  y[1:nu[j], , j] <- Directional::rsespc(nu[j], mu[j, ], theta[j, ])
  x <- y[1:nu[1], , 1]
  for (j in 2:g)  x <- rbind(x, y[1:nu[j], , j])  ## simulated data
  ## data come from the first cluster, then from the second and so on
  list(id = ina, x = x)
}
