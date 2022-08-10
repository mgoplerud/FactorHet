#ifndef eigen_utilities
#define eigen_utilities

Eigen::MatrixXd softmax_matrix(Eigen::MatrixXd X);

Eigen::SparseMatrix<double> sparse_diag(Eigen::ArrayXd m);

#endif
