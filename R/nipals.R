## NIPALS PCA
nipals_pca <- function(X, ncomp = NULL, scale = FALSE, tol = 1e-8){
  ## Some initilizations
  X <- as.matrix(X)
  nvar <- ncol(X)
  nobs <- nrow(X)

  ## Centering X matrix
  meanX <- colMeans(X)
  X <- X - rep(meanX, each = nobs)
  attr(X, "scaled:center") <- meanX
  attr(X, "name") <- deparse(match.call()$X)

  ## Setting up number of Components
  if (!is.null(ncomp)) {
    if (ncomp > min(dim(X))) stop (paste("Number of components must be less than", min(dim(X))))
  } else {
    ncomp <- min(dim(X))
  }

  W <- matrix(0, nvar, ncomp)
  T <- matrix(0, nobs, ncomp)

  ## NIPALS Loop with deflation
  for (a in 1:ncomp) {
    t <- X[, which.max(colSums(X * X))]
    t_old <- 0
    repeat{
      w <- crossprod(X, t)
      w <- w/sqrt(c(crossprod(w)))
      t <- X %*% w
      if (sum(abs(t - t_old) / t) < tol) break()
      t_old <- t
    }

    ## Deflation
    X <- X - tcrossprod(t, w)

    ## Saving
    W[, a] <- w
    T[, a] <- t
  }

  ## Naming W and T
  dimnames(W) <- list(colnames(X), paste0("Comp", 1:ncomp))
  dimnames(T) <- list(rownames(X), paste0("Comp", 1:ncomp))

  return(list(
    loadings = W,
    scores = T
  ))
}

## NIPALS Fit
nipals_pls <- function(X, Y, ncomp = NULL, scale = FALSE, tol = 1e-8) {
  ## Some initializations
  X     <- as.matrix(X)
  Y     <- as.matrix(Y)
  nobs  <- nrow(X)
  npred <- ncol(X)
  nresp <- ncol(Y)


  ## Number of components
  nc <- min(nobs, npred)
  if (!is.null(ncomp)) {
    if (ncomp > nc)
      warning(paste("Number of components can not be bigger than", min(nobs, npred)))
  } else {
    ncomp <- nc
  }

  ## Initialize some output matrices
  W <- P <- matrix(0, nrow = npred, ncol = ncomp) # Loadings and Loading Weights of Predictors
  U <- T <- matrix(0, nrow = nobs, ncol = ncomp)  # Scores of Y and X
  Q <- matrix(0, nrow = nresp, ncol = ncomp)

  ## Scaling and Centering
  meanX <- colMeans(X)
  meanY <- colMeans(Y)
  X <- X - rep(meanX, each = nobs)
  Y <- Y - rep(meanY, each = nobs)
  attr(X, "scaled:center") <- meanX
  attr(Y, "scaled:center") <- meanY
  attr(X, "name") <- deparse(match.call()$X)
  attr(Y, "name") <- deparse(match.call()$Y)

  ## Saving Original Scaled X and Y for Model Frame
  X0 <- X; Y0 <- Y;

  ## Some functions
  modelFrame <- function(which = NULL){
    if (!is.null(which)) return(switch(which, x = X0, y = Y0))
    out <- data.frame(x = I(X0), y = I(Y0))
    names(out) <- c(attr(X, "name"), attr(Y, "name"))
    return(out)
  }

  for (nc in 1:ncomp) {
    u     <- if (nresp == 1) Y else Y[, which.max(colSums(Y * Y))]
    t_old <- 0
    repeat {
      w     <- crossprod(X, u)           # X-loading weight
      w     <- w/sqrt(c(crossprod(w)))   # Normalize x-loading weight
      t     <- X %*% w                   # X-scores (unscaled)
      tt    <- t/c(crossprod(t))         # X-scores (scaled)
      q     <- crossprod(Y, tt)          # Y-loadings from scaled X-scores
      if (nresp == 1) break()
      if (sum(abs(t - t_old)/t) < tol) break()
      u     <- Y %*% q / c(crossprod(q)) # Y-Scores from scaled Y-loadings
      t_old <- t
    }
    p <- crossprod(X, tt) # X-loadings from scaled x-scores

    ## Deflation
    X <- X - tcrossprod(t, p)
    Y <- Y - tcrossprod(t, q)

    ## Saving the output
    W[, nc] <- unname(w)
    P[, nc] <- unname(p)
    Q[, nc] <- unname(q)
    T[, nc] <- unname(t)
    U[, nc] <- unname(u)
  }

  ## Naming Matrices
  colnames(W) <- colnames(P) <- colnames(T) <- colnames(U) <- colnames(Q) <- paste0("Comp", 1:ncomp)
  rownames(T) <- rownames(U) <- rownames(X)
  rownames(W) <- rownames(P) <- colnames(X)
  rownames(Q) <- colnames(Y)

  ret <- list(
    call            = match.call(),
    ncomp           = ncomp,
    nresp           = nresp,
    modelFrame      = modelFrame,
    loading.weights = W,
    loadings        = P,
    Yloadings       = Q,
    Yscores         = U,
    scores         = T
  )
  class(ret) <- "minPLS"
  return(ret)
}

## Compute Variation Explained
varX <- function(x, ncomp) UseMethod("varX", x)
varX.minPLS <- function(x, ncomp = x$ncomp){
  P <- x$loadings[, 1:ncomp]
  tsq <- sapply(1:ncomp, function(a) c(crossprod(x$scores[, a])))
  return(colSums(P * P) * tsq)
}

## Total Variation in X
totVar <- function(x) UseMethod("totVar", x)
totVar.minPLS <- function(x){
  X <- x$modelFrame("x")
  return(sum(X * X))
}

## Summary Output
summary.minPLS <- function(x, ...){
}

## Print the object
print.minPLS <- function(x, ...){
  print(x$call)
}

## Get Projection Matrix
getProjection <- function(x) UseMethod("getProjection", x)
getProjection.minPLS <- function(x, ...){
  if (x$ncomp == 1) return(x$W)
  PW <- crossprod(x$loadings, x$loading.weight)
  PWinv <- backsolve(PW, diag(x$ncomp))
  R <- x$loading.weights %*% PWinv
  colnames(R) <- paste0("Comp", 1:x$ncomp)
  return(R)
}

## Get Beta Coefficients
coef.minPLS <- function(x, ncomp = x$ncomp, ...){
  R <- getProjection(x)
  out <- lapply(1:ncomp, function(a) {
    R[, 1:a, drop = FALSE] %*% t(x$Yloadings[, 1:a, drop = FALSE])
  })
  names(out) <- paste0("Comp", 1:ncomp)
  return(out)
}

## Get Fitted Values
fitted.minPLS <- function(x, ncomp = x$ncomp, ...){
  yMeans <- attr(x$modelFrame("y"), "scaled:center")
  out <- lapply(1:ncomp, function(a){
    f <- x$scores[, 1:a, drop = F] %*% t(x$Yloadings[, 1:a, drop = F])
    f <- f + rep(yMeans, each = nrow(f))
    return(f)
  })
  names(out) <- paste0("Comp", 1:ncomp)
  return(out)
}

## Get Residuals
residuals.minPLS <- function(x, ncomp = x$ncomp, ...){
  Y <- x$modelFrame("y")
  attributes(Y) <- NULL
  out <- lapply(1:ncomp, function(a){
    Y <- Y - x$scores[, a] %*% t(x$Yloadings[, a])
    rownames(Y) <- rownames(x$scores)
    return(Y)
  })
  names(out) <- paste0("Comp", 1:ncomp)
  return(out)
}

