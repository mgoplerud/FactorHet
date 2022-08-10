#include "RcppEigen.h"
#include "eigen_utilities.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// Fast Calculation of Row-Wise SoftMax of Matrix
// [[Rcpp::export]]
Eigen::MatrixXd softmax_matrix(
    const Eigen::MatrixXd X
){
  
  Eigen::MatrixXd adj_X = X.adjoint();
  
  int N = X.rows();
  
  Eigen::MatrixXd group_Eprob(X.cols(), N);
  
  for (int i = 0; i < N; i++){
    Eigen::ArrayXd x = adj_X.col(i);
    x = x - x.maxCoeff();
    x = x.exp();
    group_Eprob.col(i) = x / x.sum();
  }
  return group_Eprob.adjoint();
}

// Sparse Diagonal
// [[Rcpp::export]]
Eigen::SparseMatrix<double> sparse_diag(
  Eigen::ArrayXd m
){
  
  int n = m.size();
  
  Eigen::SparseMatrix<double> res(n, n);
  
  std::vector<Eigen::Triplet<double> > tripletList;
  
  tripletList.reserve(n);
  
  for (int i = 0; i < n; i++) {
    
    tripletList.push_back(Eigen::Triplet<double>(i, i, m[i]));
    
  }
  
  res.setFromTriplets(tripletList.begin(), tripletList.end());
  
  return res;
}

