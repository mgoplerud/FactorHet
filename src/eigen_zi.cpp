#include "RcppEigen.h"
#include "eigen_utilities.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd logpi_adjust(
  const Eigen::Map<Eigen::MatrixXd> X,
  const Eigen::ArrayXd logpi
){
  Eigen::MatrixXd adj_X = X.adjoint();
  // Eigen::MatrixXd output(X.cols(), X.rows());
  
  for (int i = 0; i < X.rows(); i++){
    Eigen::ArrayXd xi = adj_X.col(i);
    adj_X.col(i) = xi + logpi;
  }
  return(adj_X.adjoint());
}




