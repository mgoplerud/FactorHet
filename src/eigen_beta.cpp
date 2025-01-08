#include "RcppEigen.h"
#include "eigen_utilities.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Update Beta Given Augmentation Parameters
// [[Rcpp::export]]
Eigen::MatrixXd cpp_beta(
    const int K,
    const Eigen::MappedSparseMatrix<double> X,
    const Rcpp::List E_ridge,
    const Eigen::Map<Eigen::ArrayXd> y,
    const Eigen::Map<Eigen::ArrayXd> weights,
    const Eigen::Map<Eigen::MatrixXd> E_omega,
    const Eigen::Map<Eigen::MatrixXd> obs_E_prob
){

  Eigen::MatrixXd new_beta(X.cols(), K);
  Eigen::MatrixXd adj_X = X.adjoint();
  Eigen::VectorXd adj_y = weights.array() * (y.array() - 0.5);
  
  // Loop over clusters
  for (int k = 0; k < K; k++){
    
    Eigen::VectorXd s_ik = obs_E_prob.col(k);
    
    s_ik = s_ik.cwiseProduct(adj_y);
    Eigen::SparseMatrix<double> omega_ik = sparse_diag(E_omega.col(k));
    
    Eigen::SparseMatrix<double> tau_k = E_ridge[k];
    Eigen::MatrixXd outer_X = adj_X * omega_ik * X + tau_k;
    
    new_beta.col(k) = outer_X.llt().solve(adj_X * s_ik);
  }
  
  return new_beta;
  
}

// Direct Solution Using Cpp
// [[Rcpp::export]]
Eigen::MatrixXd cpp_beta_plain(
    const int K,
    const Eigen::MappedSparseMatrix<double> X,
    const Eigen::MappedSparseMatrix<double> omega,
    const Eigen::MappedSparseMatrix<double> ridge,
    const Eigen::Map<Eigen::VectorXd> s
){
  Eigen::SparseMatrix<double> adj_X = X.adjoint();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > decomp(adj_X * omega * X + ridge);
  Eigen::MatrixXd ridge_estimates = decomp.solve(adj_X * s);
  return ridge_estimates;
}

// Conjugate Gradient Function
// [[Rcpp::export]]
List cg_custom(
    const int K,
    const Eigen::MappedSparseMatrix<double> X,
    const Eigen::Map<Eigen::MatrixXd> omega,
    const List list_ridge,
    const Eigen::Map<Eigen::VectorXd> s,
    const Eigen::Map<Eigen::MatrixXd> old_beta,
    const Eigen::Map<Eigen::MatrixXd> weights,
    const double tol,
    const int it_max = 0,
    const int low_dimension = 5
){

  bool do_weights;
  if (weights.cols() == 0){
    do_weights = false;
  }else{
    do_weights = true;
    if (weights.rows() != X.rows()){
      Rcpp::stop("weights must have either 0 cols or same dimensionality as X");
    }
  }
  
  int p_X = X.cols();
  
  Eigen::VectorXd k_it(K);
  Eigen::VectorXd k_error(K);
  Eigen::VectorXd k_convg(K);
  Eigen::MatrixXd new_beta(p_X, K);
  
  for (int k = 0; k < K; k++){

    Eigen::SparseMatrix<double> ridge = list_ridge[k];
    
    Eigen::VectorXd sqrt_omega_k = omega.col(k).cwiseSqrt();
    
    Eigen::VectorXd adj_s(X.rows());
    
    if (do_weights){
      adj_s = s.cwiseProduct(weights.col(k)).cwiseQuotient(sqrt_omega_k);
    }else{
      adj_s = s.cwiseQuotient(sqrt_omega_k);
    }
    Eigen::SparseMatrix<double> adj_X = sparse_diag(sqrt_omega_k) * X;
    
    // If low dimensional problem, just solve directly
    if (p_X <= low_dimension){
      // Eigen::MatrixXd outer_X = adj_X.adjoint() * adj_X + ridge;
      // new_beta.col(k) = outer_X.llt().solve(adj_X.adjoint() * adj_s);
      Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > decomp(adj_X.adjoint() * adj_X + ridge);
      new_beta.col(k) = decomp.solve(adj_X.adjoint() * adj_s);
      k_it(k) = 0;
      k_error(k) = 0;
      k_convg(k) = 2;
      continue;
    }
    
    Eigen::VectorXd precond_diag(p_X);
    
    for (int i = 0; i < p_X; i++){
      precond_diag(i) = 1/(adj_X.col(i).squaredNorm() + ridge.coeff(i,i));
    }
    // Eigen::SparseMatrix<double> sparse_precond = sparse_diag(precond_diag);
    
    Eigen::VectorXd beta = old_beta.col(k);

    Eigen::VectorXd residual = adj_s - adj_X * beta;
    
    double rhsNorm2 = (adj_X.adjoint() * adj_s).squaredNorm();
    Eigen::VectorXd normal_residual = adj_X.adjoint() * residual - ridge * beta;
    
    Eigen::VectorXd p = normal_residual.cwiseProduct(precond_diag);
    double absNew = normal_residual.dot(p);
    
    double threshold = tol * tol * rhsNorm2;
    double tol_error = 0;
    
    int internal_cg_it;
    if (it_max == 0){
      internal_cg_it = p_X;
    }else{
      if (p_X < it_max){
        internal_cg_it = p_X;
      }else{
        internal_cg_it = it_max;
      }
    }
    
    // Do conjugate gradient with normal equation
    // (X^T O X + R) beta = X^T y but never form X^T O X
    int it = 0;
    bool convg = false;
    
    while (it < internal_cg_it){
      
      Eigen::VectorXd tmp = adj_X * p;
      
      double alpha_denom = tmp.squaredNorm() + p.adjoint() * ridge * p;
      double alpha = absNew / alpha_denom;
      
      beta += alpha * p;
      residual -= alpha * tmp;
      
      Eigen::VectorXd normal_residual = adj_X.adjoint() * residual - ridge * beta;
      
      double residualNorm2 = normal_residual.squaredNorm();
      tol_error = std::sqrt(residualNorm2 / rhsNorm2);
      
      if (residualNorm2 < threshold){
        convg = true;
        break;
      }
      
      Eigen::VectorXd z = normal_residual.cwiseProduct(precond_diag);
      double absOld = absNew;
      absNew = normal_residual.dot(z);
      double step_beta = absNew / absOld;
      
      p = z + step_beta * p;
      
      it++;
    }
    
    new_beta.col(k) = beta;
    k_it(k) = it;
    k_error(k) = tol_error;
    k_convg(k) = convg;
    
  }
  
  
  
  return List::create(
    Rcpp::Named("beta") = new_beta,
    Rcpp::Named("iter") = k_it,
    Rcpp::Named("error") = k_error,
    Rcpp::Named("converged") = k_convg
  );
  
}