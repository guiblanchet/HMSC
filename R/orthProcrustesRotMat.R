##' Orthogonal procrustean rotation matrix
##' 
##' @param X input matrix
##' @param Z target matrix
##' @return the rotation matrix
##' @export
orthProcrustesRotMat <- function(X, Z) {
    sol <- svd(crossprod(Z, X))
    return(sol$v %*% t(sol$u))
}
