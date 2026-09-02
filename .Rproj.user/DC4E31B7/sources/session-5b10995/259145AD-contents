mixkent.mle <- function(x, g = 2, tol = 1e-4) {

  fun2 <- function(wlika, rswlika, G, param, x, g, lika) {

    wij <- wlika / rswlika  ## weights
    pj <- Rfast::colmeans(wij)   # PANOS

    for (j in 1:g) {
      mod <- .wkent(x, w = wij[, j] )
      G[[ j ]] <- mod$G
      param[j, ] <- mod$param
      lika[, j] <- Directional::dkent(x, G[[ j ]], param[j, ], logden = TRUE) + log( pj[j] )
    }

    wlika <-  exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika/rswlika

    list(wij = wij, G = G, param = param, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  n <- dim(x)[1]
  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  G <- list()
  param <- matrix(nrow = g, ncol = 3)

  runtime <- proc.time()
  ## Step 1

  if ( g > 1 ) {
    #cl <- kmeans(x, g, nstart = n.start)$cl
    y <- Directional::euclid.inv(x)
    mod <- mixture::gpcm(y, G = g, mnames = NULL, start = 2, mmax = 100)
    cl <- Rfast::rowMaxs(mod$z)
  } else  cl <- rep(1, n)
  wij <- tabulate(cl)

  while ( min(wij) <= 10 ) {
    g <- g - 1
    lika <- matrix(nrow = n, ncol = g)
    param <- matrix(nrow = g, ncol = 3)
    #cl <- kmeans(x, g, nstart = n.start)$cl
    mod <- mixture::gpcm(y, G = g, mnames = NULL, start = 2, mmax = 100)
    cl <- Rfast::rowMaxs(mod$z)
    wij <- tabulate(cl)
  }

  for ( j in 1:g ) {
    mod <- Directional::kent.mle(x[cl == j, ])
    G[[ j ]] <- mod$G
    param[j, ] <- mod$param
    lika[, j] <-  Directional::dkent(x, G[[ j ]], param[j, ])
  }

  wlika <- exp(lika)
  rswlika <- Rfast::rowsums(wlika)

  ep <- fun2(wlika, rswlika, G, param, x, g, lika)
  lik[1] <- ep$lik
  ep2 <- fun2(ep$wlika, ep$rswlika, ep$G, ep$param, x, g, ep$lika)
  lik[2] <- ep2$lik

  i <- 2
  while ( lik[i] - lik[i - 1] > tol ) {
    i <- i + 1
    ep <- ep2
    ep2 <- fun2(ep$wlika, ep$rswlika, ep$G, ep$param, x, g, ep$lika)
    lik[i] <- ep2$lik
  }
  res <- ep2
  if ( ep$lik > ep2$lik )  res <- ep

  pj <- Rfast::colmeans(res$wij)
  loglik <- res$lik   #sum( log( Rfast::colsums( pj * t( exp( res$lika ) ) ) ) )
  ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
  param <- cbind(pj, res$param)
  G <- res$G
  runtime <- proc.time() - runtime
  colnames(param) <- c( "probs", "kappa", "beta", "psi" )
  rownames(param) <- paste("Cluster", 1:g, sep = " ")

  list( param = param, G = G, loglik = loglik,
        probs = res$wij, pred = ta, iters = i, runtime = runtime )

}




.wkent <- function(x, w) {
 
  W <- sum(w)            
  xbar <- Rfast::eachcol.apply(x, w)
  xbar <- xbar / sqrt( sum(xbar^2) )
  S <- crossprod(x, weights * x) / W   
  u <- c( acos(xbar[1]), ( atan(xbar[3] / xbar[2]) + pi * I(xbar[2] < 0) ) %% (2 * pi) )
  theta <- u[1]   ;   phi <- u[2]
  costheta <- cos(theta) ;   sintheta <- sin(theta)
  cosphi <- cos(phi)     ;   sinphi <- sin(phi)
  H <- matrix( c(costheta, sintheta * cosphi, sintheta * sinphi,
                 -sintheta, costheta * cosphi, costheta * sinphi,
                 0, -sinphi, cosphi), ncol = 3)
  B <- crossprod(H, S) %*% H
  psi <- 0.5 * atan(2 * B[2, 3] / (B[2, 2] - B[3, 3]))
  K <- matrix( c(1, 0, 0,
                 0, cos(psi), sin(psi),
                 0, -sin(psi), cos(psi)), ncol = 3 )
  G <- H %*% K
  colnames(G) <- c("mean", "major", "minor")
  
  xg <- x %*% G
  xg1 <- sum(w * xg[, 1])
  xg2 <- sum(w * xg[, 2]^2)
  xg3 <- sum(w * xg[, 3]^2)
  
    mle <- function(para) {
      k <- para[1] ; b <- para[2]
      ckb <- Directional::kent.logcon(k, b)   
      W * ckb - k * xg1 - b * (xg2 - xg3)     
    }
  
  ini_kappa <- Rfast::vmf.mle(x)$kappa
  ini <- c(ini_kappa, ini_kappa / 2.1) 
  qa <- optim(ini, mle, control = list(maxit = 5000))
  para <- qa$par
  k <- para[1]   ;   b <- para[2]
  ckb <- Directional::kent.logcon(k, b)
  lik <-  - W * ckb + k * xg1 + b * (xg2 - xg3)
  param <- c(k, b, psi)
  names(param) <- c("kappa", "beta", "psi") 
  list(G = G, param = param, logcon = ckb, loglik = lik)
}
