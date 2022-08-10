#include "RcppEigen.h"
#include "eigen_utilities.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// C++ implementation of tr([X^T O X + R]^{-1} (X^T O X)). 
// For logistic regression
// [[Rcpp::export]]
double trace_df_cpp(
  const Eigen::MappedSparseMatrix<double> X,
  const Eigen::Map<Eigen::MatrixXd> beta,
  const Eigen::MappedSparseMatrix<double> ridge
){
  
  int p = X.cols();
  double df_kc;


  // Get p = 1/(1 + exp(-x^T beta)) = exp(x^T beta) / (1 + exp(x^T beta))
  Eigen::ArrayXd xb = (X * beta).array();
  Eigen::ArrayXd pr = 1/(1 + (-xb).exp());
  // Get the Hessian: X^T Diag[p * (1-p)] * X
  Eigen::ArrayXd meat = pr - pr * pr;
  
  Eigen::SparseMatrix<double> Hess = X.adjoint() * sparse_diag(meat) * X;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(Hess + ridge);
  
  //As (H + R)^{-1} H may be very large, avoid every explicitly computing it
  //as we only need the *diagonal* elements Thus, perform a solve for each
  //column of the Hessian and extract only the corresponding diagonal element
  Eigen::VectorXd diag_j(p);
  for (int j = 0; j < p; j++){
    Eigen::VectorXd inv_j = llt.solve(Hess.col(j));
    diag_j[j] = inv_j(j);
  }
  
  df_kc = diag_j.sum();

  return df_kc;
  
}

// [[Rcpp::export]]
double trace_AinvB_sparse(
    const Eigen::MappedSparseMatrix<double> A,
    const Eigen::MappedSparseMatrix<double> B
){
  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(A);
  
  double out_trace;
  int p = B.rows();
  
  Eigen::VectorXd diag_out(p);
  for (int j = 0; j < p; j++){
    Eigen::VectorXd inv_j = llt.solve(B.col(j));
    diag_out[j] = inv_j(j);
  }
  
  out_trace = diag_out.sum();
  
  return out_trace;
}

// [[Rcpp::export]]
double trace_AinvB(
    const Eigen::MappedSparseMatrix<double> A,
    const Eigen::Map<Eigen::MatrixXd> B
){
  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(A);
  
  double out_trace;
  int p = B.rows();
  
  Eigen::VectorXd diag_out(p);
  for (int j = 0; j < p; j++){
    Eigen::VectorXd inv_j = llt.solve(B.col(j));
    diag_out[j] = inv_j(j);
  }
  out_trace = diag_out.sum();

  return out_trace;
}

// [[Rcpp::export]]
double fast_sparse_df(
    const Eigen::MappedSparseMatrix<double> LL,
    const Eigen::MappedSparseMatrix<double> PR
){
  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > llt;
  llt.compute(LL + PR);
  
  double out_trace;
  int p = LL.rows();
  
  Eigen::VectorXd diag_out(p);
  for (int j = 0; j < p; j++){
    Eigen::VectorXd inv_j = llt.solve(LL.col(j));
    diag_out[j] = inv_j(j);
  }
  
  out_trace = diag_out.sum();
  
  return out_trace;
}
