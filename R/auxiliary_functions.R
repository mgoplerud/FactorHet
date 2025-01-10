#' Rank of Matrix
#' 
#' Calculate rank of (sparse) matrix by a QR decomposition.
#' @keywords internal
#' @param X Sparse Matrix with size N by p.
#' @param outer Calculate rank of X^T X. Often much faster if N >> p.
#' @import Matrix
rank_via_null <- function(X, outer=FALSE){
  dim_X <- dim(X)  
  if (outer){
    X <- crossprod(X)
    X <- as(X, 'dgCMatrix')
  }else{
    X <- as(X, 'dgCMatrix')
  }
  rank_X <- rank_sparse(X)
  return(rank_X)
}

# A simple alternative for "print"
#' @importFrom utils capture.output
easy_message <- function(x){
  message(paste0(capture.output(x), collapse = "\n"))
}