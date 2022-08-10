#include "RcppEigen.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Fast Calculation of Basis for Null Space.
// [[Rcpp::export]]
Eigen::SparseMatrix<double> calculate_nullspace_basis(
    const Eigen::MappedSparseMatrix<double> X
){
  //Perform QR
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > obj;
  //obj.setPivotThreshold(tol);
  obj.compute(X);
  int rank = obj.rank();
  int ncolX = X.cols();
  
  // Address two edge cases
  if (rank == ncolX){
    Eigen::SparseMatrix<double> empty(ncolX,0);
    //Eigen::SparseMatrix<double> empty(ncolX, 1);
    return empty;
  }
  if (rank == 0){
    Eigen::SparseMatrix<double> I(ncolX, ncolX);
    I.setIdentity();
    return I;
  }
  
  Eigen::SparseMatrix<double,Eigen::RowMajor> Rr = obj.matrixR();  // row-major, sorted
  Eigen::SparseMatrix<double> Rc = Rr;            // column-major, sorted
  Eigen::SparseMatrix<double> R1 = Rc.topLeftCorner(rank, rank);
  Eigen::SparseMatrix<double> R2 = Rc.topRightCorner(rank, ncolX - rank);
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >::PermutationType Pmat = obj.colsPermutation();
  
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solve_square;
  solve_square.compute(R1);
  
  Eigen::SparseMatrix<double> upper = solve_square.solve(R2);
  
  //Create a matrix to flip the sign of R1^{-1} R2.
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(upper.rows());
  for (int i = 0; i < upper.rows(); ++i){
    tripletList.push_back(T(i,i,-1));
  }
  Eigen::SparseMatrix<double> negative_identity(upper.rows(), upper.rows());
  negative_identity.setFromTriplets(tripletList.begin(), tripletList.end());
  upper = negative_identity * upper;
  
  int size_i = ncolX - rank;
  Eigen::SparseMatrix<double> I(size_i, size_i);
  I.setIdentity();
  
  Eigen::SparseMatrix<double> out(ncolX - rank, ncolX);
  out.leftCols(rank) = upper.adjoint();
  out.rightCols(size_i) = I;
  Eigen::SparseMatrix<double> basis = Pmat * out.adjoint();
  return basis;
}


// Rank of matrix from QR
// [[Rcpp::export]]
int rank_sparse(
    const Eigen::MappedSparseMatrix<double> X
){
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > obj;
  //obj.setPivotThreshold(tol);
  obj.compute(X);
  int rank = obj.rank();
  return rank;
}