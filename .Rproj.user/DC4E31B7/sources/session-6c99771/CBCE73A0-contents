mixesag.mle <- function(x, g = 2, n.start = 10, tol = 1e-4) {

  fun <- function(param, z, w, I3) {
    m <- param[1:3]
    gam1 <- param[4]   ;   gam2 <- param[5]
    heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
    m0 <- sqrt( m[2]^2 + m[3]^2 )
    rl <- sum(m^2)
    x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / m0 / sqrt(rl)
    x2b <- c( 0, -m[3], m[2] )/m0
    T1 <- tcrossprod( x1b )   ;    T2 <- tcrossprod( x2b )
    T12 <- tcrossprod( x1b, x2b )
    vinv <- I3 + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + heta * ( T1 + T2 )
    g2 <- colSums( m * z )
    g1 <- colSums( z * crossprod(vinv, z) )
    a <- g2 / sqrt(g1)
    a2 <- a^2
    M2 <- ( 1 + a2 ) * pnorm(a) + a * dnorm(a)
    - 0.5 * sum(w * a2) + 0.5 * sum(w) * rl + 1.5 * sum( w * log(g1) ) - sum( w * log(M2) )
  }

  fun2 <- function(wlika, rswlika, param, z, I3, g, lika) {

    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS

    for (j in 1:g) {
      a1 <- optim(param[j, ], fun, z = z, w = wij[, j], I3 = I3, control = list(maxit = 5000), )
      a2 <- optim(a1$par, fun, z = z, w = wij[, j], I3 = I3, control = list(maxit = 5000), method = "BFGS")
      while (a1$value - a2$value > 1e-05) {
        a1 <- a2
        a2 <- optim(a1$par, fun, z = z, w = wij[, j], I3 = I3, control = list(maxit = 5000) )
      }
      param[j, ] <- a2$par
      m <- param[j, 1:3]
      gam1 <- param[j, 4]   ;   gam2 <- param[j, 5]
      heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
      m0 <- sqrt( m[2]^2 + m[3]^2 )
      rl <- sum(m^2)
      x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / m0 / sqrt(rl)
      x2b <- c( 0, -m[3], m[2] )/m0
      T1 <- tcrossprod( x1b )   ;   T2 <- tcrossprod( x2b )
      T12 <- tcrossprod( x1b, x2b )
      vinv <- I3 + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + heta * ( T1 + T2 )
      g2 <- colSums( m * z )
      g1 <- colSums( z * crossprod(vinv, z) )
      a <- g2 / sqrt(g1)
      a2 <- a^2
      M2 <- ( 1 + a2 ) * pnorm(a) + a * dnorm(a)
      lika[, j] <-  - ( - 0.5 * a2 + 0.5 * rl + 1.5 * log(g1) - log(M2) ) + log( pj[j] )
    }

    wlika <-  exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika/rswlika

    list(wij = wij, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  n <- dim(x)[1]
  I3 <- diag(3)
  z <- t(x)

  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  param <- matrix(nrow = g, ncol = 5)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    cl <- kmeans(x, g, nstart = n.start)$cl
  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 10 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 5)
    cl <- kmeans(x, g, nstart = n.start)$cl
    wij <- tabulate(cl)
  }

  for ( j in 1:g ) {
    mod <- Directional::esag.mle(x[cl == j, ])
    param[j, ] <- c(mod$mu, mod$gam)
    m <- param[j, 1:3]
    gam1 <- param[j, 4]   ;   gam2 <- param[j, 5]
    heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
    m0 <- sqrt( m[2]^2 + m[3]^2 )
    rl <- sum(m^2)
    x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / m0 / sqrt(rl)
    x2b <- c( 0, -m[3], m[2] )/m0
    T1 <- tcrossprod( x1b )   ;   T2 <- tcrossprod( x2b )
    T12 <- tcrossprod( x1b, x2b )
    vinv <- I3 + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + heta * ( T1 + T2 )
    g2 <- colSums( m * z )
    g1 <- colSums( z * crossprod(vinv, z) )
    a <- g2 / sqrt(g1)
    a2 <- a^2
    M2 <- ( 1 + a2 ) * pnorm(a) + a * dnorm(a)
    lika[, j] <-  - ( - 0.5 * a2 + 0.5 * rl + 1.5 * log(g1) - log(M2) )
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  ep <- fun2(wlika, rswlika, param, z, I3, g, lika)
  lik[1] <- ep$lik
  ep2 <- fun2(ep$wlika, ep$rswlika, ep$param, z, I3, g, ep$lika)
  lik[2] <- ep2$lik

  i <- 2
  while ( lik[i] - lik[i - 1] > tol ) {
    i <- i + 1
    ep <- ep2
    ep2 <- fun2(ep$wlika, ep$rswlika, ep$param, z, I3, g, ep$lika)
    lik[i] <- ep2$lik
  }
  res <- ep2
  if ( ep$lik > ep2$lik )  res <- ep

  pj <- Rfast::colmeans(res$wij)
  loglik <- res$lik   #sum( log( Rfast::colsums( pj * t( exp( res$lika ) ) ) ) )
  ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
  param <- cbind(pj, res$param)
  runtime <- proc.time() - runtime
  colnames(param) <- c( "probs", paste("mu", 1:3, sep = ""), "gamma1", "gamma2" )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")

  list( param = param, loglik = loglik - n * log(2 * pi),
        probs = res$wij, pred = ta, iters = i, runtime = runtime )

}
