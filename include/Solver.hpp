#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF

#include "CSR_Matrix.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>

class Solver
{
	public:
		int n;
		int p;
		CSRMatrix A;
		std::vector<double> b;
		std::vector<double> u;
		Eigen::SparseMatrix<double> EigenA;
		Eigen::VectorXd Eigenb;

		// constructor
		Solver(int n, int p,CSRMatrix A, std::vector<double> b, std::vector<double> u);

		// destructor
		~Solver();

		// setup system
		void setupEigen();
		// solve system
		std::vector<double> solveEigen();
};

#endif