# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

minPLS <- function(X, Y, ncomp, ...) {
  ## Collecting Information
  Y <- as.matrix(Y)
  nobs <- nrow(X)
  npred <- ncol(X)
  nresp <- ncol(Y)

  ## Scale X and Y
  xmeans <- apply(X, 2, mean)
  ymeans <- apply(Y, 2, mean)
  X <- X - rep(xmeans, each = nobs)
  Y <- Y - rep(ymeans, each = nobs)

  ## Take a vector U from Y
  U <- apply(Y, 1, max)


}
