#include "Solver.hpp"
#include <iostream>

// constructor
Solver::Solver(int n, int p, CSRMatrix A, std::vector<double> b, std::vector<double> u)
 : n(n), p(p), A(A), b(b), u(u)
{
}

// destructor
Solver::~Solver()
{
}

void Solver::setupEigen()
{
	int noNodes = b.size();
	EigenA.resize(noNodes, noNodes);
	Eigenb.resize(noNodes);

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i=0; i<noNodes; i++)
	{
		for (int j=A.row_start[i]; j<A.row_start[i+1]; j++)
		{
			triplets.emplace_back(Eigen::Triplet<double>(i, A.col_no[j], A.entries[j]));
		}
	}

	EigenA.setFromTriplets(triplets.begin(), triplets.end());
	EigenA.makeCompressed();

	for (int i=0; i<noNodes; i++)
	{
		Eigenb(i) = b[i];
	}
}

std::vector<double> Solver::solveEigen()
{
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(EigenA);

	if (solver.info() != Eigen::Success)
	{
		throw std::runtime_error("Eigen factorisation failed");
	}

	Eigen::VectorXd solution = solver.solve(Eigenb);

	for (int i=0; i<static_cast<int>(b.size()); i++)
	{
		u[i] = solution(i);
	}

	return u;
}