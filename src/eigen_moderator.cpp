#include "RcppEigen.h"
#include "eigen_utilities.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;


// Moderator Hessian and Gradient
// [[Rcpp::export]]
Rcpp::List cpp_moderator_deriv(
    const Eigen::Map<Eigen::MatrixXd> W,
    const Eigen::Map<Eigen::MatrixXd> phi,
    const bool do_grad = true,
    const bool do_hess = true,
    const bool do_ln_hess = true,
    const bool moderator_grad = true
){
  
  int K = phi.rows();
  int size_W = W.cols();
  int nrow_W = W.rows();
  Eigen::MatrixXd adj_W = W.adjoint();
  Eigen::MatrixXd pi = W * phi.adjoint();
  pi = softmax_matrix(pi);

  // Eigen::ArrayXd pi_bar = ((norm_weights.asDiagonal()) * pi_ik).colwise().sum().array();
  
  Eigen::ArrayXd pi_bar = pi.colwise().mean().array();
  
  Eigen::VectorXd neg_sq_pibar = pi_bar.pow(-2);
  
  Rcpp::List grad_output(K);
  Rcpp::List hessian_output(K);
  Rcpp::List hessian_ln_output(K);

  for (int kprime = 0; kprime < K; kprime++){
    Eigen::MatrixXd grad_mat_kprime(size_W, K);
    Eigen::VectorXd pi_kprime = pi.col(kprime);
    Rcpp::List k_inner(K - 1);
    Rcpp::List k_ln_inner(K - 1);
    if (do_grad or do_ln_hess or moderator_grad){
      for (int k = 1; k < K; k++){
        // Rcout << kprime << "-" << k << std::endl;
        Eigen::VectorXd pi_k = pi.col(k);
        Eigen::VectorXd grad_mat_kkprime(size_W);
        Eigen::VectorXd grad_matkkprime_weight(nrow_W);
        if (k == kprime){
          grad_matkkprime_weight.array() = pi_k.array() - pi_k.array() * pi_kprime.array();
        }else{
          grad_matkkprime_weight.array() = -pi_k.array() * pi_kprime.array();
        }
        // Rcout << grad_matkkprime_weight.mean() << std::endl; 
        grad_mat_kkprime = (grad_matkkprime_weight.asDiagonal() * W).colwise().mean();
        grad_mat_kprime.col(k) = grad_mat_kkprime;
      }
    }
    
    for (int k =1; k < K; k++){
      Eigen::VectorXd pi_k = pi.col(k);
      Eigen::VectorXd out_weight = - pi_k;
      out_weight.array() += (kprime == k);
      
      Rcpp::List m_inner(K - 1);
      Rcpp::List m_ln_inner(K - 1);
      for (int m = 1; m < K; m++){
        Eigen::VectorXd pi_m = pi.col(m);
        Eigen::MatrixXd hess_out_km(size_W, size_W);
        if (do_hess or do_ln_hess or moderator_grad){
          Eigen::VectorXd hess_weight_km = - pi_m;
          hess_weight_km.array() += (kprime == m);
          hess_weight_km = out_weight * hess_weight_km;          
          
          Eigen::VectorXd hess_weight_adjust = -pi_m;
          hess_weight_adjust.array() += (m == k);
          
          hess_weight_km.array() -= hess_weight_adjust.array() * pi_k.array();
          hess_weight_km.array() = hess_weight_km.array() * pi_kprime.array()/nrow_W;
          
          hess_out_km = (adj_W * (hess_weight_km.asDiagonal()) ) * W;
          m_inner(m - 1) = hess_out_km;
        }
        if (do_ln_hess or moderator_grad){
          Eigen::MatrixXd ln_hess_out_km = hess_out_km * 1/pi_bar[kprime];
          ln_hess_out_km -= neg_sq_pibar[kprime] * (grad_mat_kprime.col(k) * grad_mat_kprime.col(m).transpose());
          m_ln_inner(m-1) = ln_hess_out_km;
        }
      }
      k_ln_inner(k - 1) = m_ln_inner;  
      k_inner(k - 1) = m_inner;
    }
    hessian_ln_output(kprime) = k_ln_inner;
    hessian_output(kprime) = k_inner;
    grad_output(kprime) = grad_mat_kprime;
  }
  
  return List::create(
    Rcpp::Named("grad_pibar") = grad_output,
    Rcpp::Named("hessian_pibar") = hessian_output,
    Rcpp::Named("hessian_lnpibar") = hessian_ln_output
  );
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_gradient_phi(
    const Eigen::Map<Eigen::VectorXd> par, 
    const int K, 
    const Eigen::Map<Eigen::VectorXd> norm_weights, 
    const Eigen::Map<Eigen::VectorXd> weights_W, 
    const Eigen::Map<Eigen::MatrixXd> group_E_prob, 
    const Eigen::Map<Eigen::MatrixXd> W, 
    const Eigen::Map<Eigen::MatrixXd> ridge_penalty, 
    const double gamma, 
    const double rank_F, 
    const double power_pi, 
    const Eigen::ArrayXd b_r, 
    const double lambda,
    Nullable<NumericVector> sampling_weights = R_NilValue
){
  int size_W = W.cols();
  int n_W = W.rows();
  
  Eigen::Map<const Eigen::MatrixXd> phi(par.data(), K - 1, size_W);
  
  Eigen::MatrixXd pi_ik = Eigen::MatrixXd::Zero(n_W, K);
  pi_ik.rightCols(K-1) = W * phi.adjoint();
  pi_ik = softmax_matrix(pi_ik);
  
  Eigen::VectorXd pi_bar = ((norm_weights.asDiagonal()) * pi_ik).colwise().sum();
  
  Eigen::MatrixXd grad_ll_phi(K - 1, size_W);
  Eigen::MatrixXd grad_reg_phi(K - 1, size_W);
  Eigen::MatrixXd output(K, size_W);

  for (int k = 1; k < K; k++){
    Eigen::VectorXd weight_ll(n_W);
    
    weight_ll = group_E_prob.col(k) - pi_ik.col(k);
    weight_ll.array() = weights_W.array() * weight_ll.array();

    grad_ll_phi.row(k-1) = ( weight_ll.asDiagonal() * W).colwise().sum();
  }
  
  for (int i = 0; i < size_W; i++){
    grad_reg_phi.col(i) = -ridge_penalty * phi.col(i);
  }

  if (gamma > 0){
    Eigen::ArrayXd weight_k(K);
    
    if (lambda > 0){
      weight_k = pi_bar.array().pow(-1) * gamma * rank_F - 
        lambda * gamma * pi_bar.array().pow(gamma - 1) * b_r;
    }else{
      weight_k = pi_bar.array().pow(-1) * gamma * rank_F - 
        gamma * pi_bar.array().pow(gamma - 1) * b_r;
    }
    
    Eigen::MatrixXd array_weights(K, size_W);
    // Rcpp::List array_list(K);
    
    for (int i = 1; i < K; i++){
      
      Eigen::VectorXd pi_i = pi_ik.col(i);
      Eigen::MatrixXd al_k(K, size_W);
      
      for (int j = 0; j < K; j++){
        
        Eigen::VectorXd pi_j = pi_ik.col(j);
        Eigen::VectorXd cross_ij(n_W);
        if (i == j){
          cross_ij.array() = pi_i.array() * (1- pi_j.array());
        }else{
          cross_ij.array() = -pi_i.array() * pi_j.array();
        }
        
        cross_ij.array() = cross_ij.array() * weights_W.array();
        
        Eigen::VectorXd weight_cross = (cross_ij.asDiagonal() * W).colwise().mean();
        weight_cross.array() = weight_k(j)  * weight_cross.array();
        al_k.row(j) = weight_cross;
      }
      array_weights.row(i) = al_k.colwise().sum();
      // array_list(i) = al_k;
    }
    output = grad_ll_phi + grad_reg_phi + array_weights.bottomRows(K-1);
  }else{
    output = grad_ll_phi + grad_reg_phi;
  }
  Eigen::Map<Eigen::VectorXd> flat_output(output.data(), output.size());

  return flat_output;
}

// [[Rcpp::export]]
double cpp_obj_phi(
    const Eigen::Map<Eigen::VectorXd> par, 
    const int K, 
    const Eigen::Map<Eigen::VectorXd> norm_weights, 
    const Eigen::Map<Eigen::VectorXd> weights_W, 
    const Eigen::Map<Eigen::MatrixXd> W, 
    const Eigen::Map<Eigen::MatrixXd> group_E_prob, 
    const Eigen::Map<Eigen::MatrixXd> ridge_penalty, 
    const double gamma, 
    const double rank_F, 
    const double power_pi, 
    const Eigen::ArrayXd b_r, 
    const double lambda
){
  // Nullable<Eigen::VectorXd> sampling_weights = R_NilValue, 
  int size_W = W.cols();
  int n_W = W.rows();
  
  Eigen::Map<const Eigen::MatrixXd> phi(par.data(), K - 1, size_W);
  
  Eigen::MatrixXd pi_ik = Eigen::MatrixXd::Zero(n_W, K);
  pi_ik.rightCols(K-1) = W * phi.adjoint();
  pi_ik = softmax_matrix(pi_ik);
  
  Eigen::ArrayXd pi_bar = ((norm_weights.asDiagonal()) * pi_ik).colwise().sum().array();

  double loglik_zik = ((weights_W.asDiagonal() * group_E_prob).array() * pi_ik.array().log()).sum();

  double regularization_phi = -0.5 * (phi.adjoint() * ridge_penalty * phi).diagonal().sum();
  
  double regularization_beta = 0;
  regularization_beta += pi_bar.log().sum() * gamma * rank_F;
  
  if (lambda > 0){
    regularization_beta -= lambda * (pi_bar.pow(power_pi) * b_r.array()).sum();
  }else{
    regularization_beta -= (pi_bar.pow(power_pi) * b_r.array()).sum();
  }
  
  double objective = loglik_zik + regularization_phi + regularization_beta;
  return objective;      
}