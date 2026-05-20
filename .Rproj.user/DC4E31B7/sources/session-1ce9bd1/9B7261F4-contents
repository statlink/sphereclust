# Geodesic maps (Fletcher 2004 style) - which doesn't do sphere-to-sphere
# Subsphere projection (Jung-Dryden-Marron 2012 style)

# Proper PGA using exponential/logarithmic maps on sphere
geodesic_pca_sphere <- function(x, k = 2) {
  # X: n × d on S^{d-1}
  # Returns projection onto S^k via geodesic PCA
  dm <- dim(x)  ;  n <- dm[1]  ;  d <- dm[2]
  # Step 1: Compute Fréchet mean (intrinsic mean on sphere)
  mu <- .frechet_mean_sphere(x)
  # Step 2: Logarithmic map: sphere → tangent space at mu
  # Log map: log_mu(x) = (angle between mu and x) * (x - <x,mu>mu) / ||x - <x,mu>mu||
  V <- .log_map_sphere(x, mu)  # n × d tangent vectors
  # Step 3: PCA in tangent space
  svd_res <- svd(V, nu = 0, nv = k)
  # Top k principal directions in tangent space
  U_k <- svd_res$v[, 1:k]  # d × k
  # Step 4: Project tangent vectors onto k-dimensional subspace

  # V_proj <- V %*% U_k %*% t(U_k)  # n × d (but rank k)

  # Step 5: Exponential map back to sphere
  # But exp_mu(V_proj) gives points on S^{d-1}, not S^k
  # To get S^k, we need to embed via the k principal geodesic directions
  # Alternative: coordinates in geodesic subspace

  #coords <- V %*% U_k  # n × k coordinates

  # Build k+1 dimensional embedding
  # First direction: mean mu
  # Next k directions: principal geodesics
  # Orthonormalize [mu, U_k] to get (k+1) directions
  basis <- qr.Q( qr( cbind(mu, U_k) ) )  # d × (k+1)
  # Project data onto this (k+1)-dim subspace and normalize
  x_proj <- x %*% basis
  x_proj / sqrt( Rfast::rowsums(x_proj^2) )
  #list(x_reduced = x_sphere, basis = basis,
  #eigenvalues = svd_res$d[1:k]^2, mean = mu, tangent_coords = coords)
}


.frechet_mean_sphere <- function(x, tol = 1e-8, maxiter = 100) {
  mu <- Rfast::colmeans(x)
  mu <- mu / sqrt( sum(mu^2))

  for ( iter in 1:maxiter ) {
    # Tangent vectors from mu to each point
    V <- .log_map_sphere(x, mu)
    # Mean tangent vector
    v_mean <- Rfast::colmeans(V)
    # Step size
    step_norm <- sqrt( sum(v_mean^2) )
    if ( step_norm < tol ) break
    # Exponential map: move mu along v_mean
    mu <- exp_map_sphere(v_mean, mu)
  }

  mu
}


# Logarithmic map: S^{d-1} → T_mu S^{d-1}
.log_map_sphere <- function(x, mu) {
  # x: n × d points on sphere
  # mu: d-vector, base point on sphere
  n <- dim(x)[1]
  # Inner products
  dots <- c( x %*% mu )
  dots <- pmax (pmin(dots, 1), -1 )  # clip to [-1,1]
  # Geodesic distance (angle)
  theta <- acos(dots)
  # Direction: (x - <x,mu>mu) normalized
  x_ortho <- x - outer(dots, mu)
  norms <- sqrt( Rfast::rowsums(x_ortho^2) )
  # Handle points at mu (zero distance)
  norms[norms < 1e-10] <- 1
  theta[theta < 1e-10] <- 0
  # Log map: theta * direction
  V <- (theta / norms) * x_ortho
  V[theta < 1e-10, ] <- 0  # zero vector at mu
  V
}

# Exponential map: T_mu S^{d-1} → S^{d-1}
.exp_map_sphere <- function(v, mu) {
  # v: tangent vector at mu
  # mu: base point on sphere
  v_norm <- sqrt( sum(v^2) )
  if ( v_norm < 1e-10 ) {
    return(mu)
  }
  # Geodesic: cos(t)*mu + sin(t)*(v/||v||)
  # At t = ||v||:
  x <- cos(v_norm) * mu + sin(v_norm) * (v / v_norm)
  x / sqrt( sum(x^2) )  # renormalize
}


