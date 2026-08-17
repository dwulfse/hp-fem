#include "FE_Solution.hpp"
#include "FE_Mesh1D.hpp"
#include "FE_Mesh2D.hpp"
#include "GaussQuadrature.hpp"
#include "GaussQuadrature2D.hpp"
#include "Element1D.hpp"
#include "Element2D.hpp"
#include "Helper.hpp"
#include "Solver.hpp"

#include <iostream>

// constructor
FE_Solution::FE_Solution(int n, int p, int d)
 : n(n), p(p), d(d), globalStiffness()
{
	// decide 1D or 2D mesh based on dimension
	if (d == 1)
	{
		mesh = std::make_unique<FE_Mesh1D>(n, p, d);
		globalLoad.resize(p*n+1);
		solution.resize(p*n+1);
	}
	else if (d == 2)
	{
		mesh = std::make_unique<FE_Mesh2D>(n, p, d);
	}
}

// destructor
FE_Solution::~FE_Solution()
{
}

// controller function for whole program
std::vector<double> FE_Solution::solve(double (*f)(const std::vector<double>&),
									   std::string fileNameNoExt,
									   double u0, double u1,
									   bool boundary_u0, bool boundary_u1)
{
	// mesh initialied in constructor
	// construct mesh from file if 2D, else construct based on n and p
	mesh->constructMesh(d == 2 ? fileNameNoExt : "");
	// allocate space for CSR stiffness matrix based on element DoFs
	mesh->allocateStiffness();
	// construct global stiffness matrix from local element stiffness matrices
	mesh->assembleStiffnessMatrix();
	// construct global load vector from local element load vectors
	mesh->assembleLoadVector(f);
	// apply boundary conditions,
	// set left and right boundaries in 1D,
	// or edge boundaries in 2D
	mesh->applyBoundaryConditions(u0, u1, boundary_u0, boundary_u1);

	// resize load and solution vectors based on number of DoFs
	int nDoF = mesh->getNoNodes();
	globalLoad.resize(nDoF);
	solution.resize(nDoF);

	globalStiffness = mesh->stiffness;
	globalLoad = mesh->load;
	n = mesh->n;

	// create Solver object and solve system
	Solver solver(n, p, globalStiffness, globalLoad, solution);

	solver.setupEigen();
	solution = solver.solveEigen();

	return solution;
}

std::vector<double> FE_Solution::solveSemilinear(double (*f)(const std::vector<double> &),
												 int q, double alpha, double lambda,
												 int maxIter, double tol,
												 std::string fileNameNoExt,
												 double u0, double u1,
												 bool boundary_u0, bool boundary_u1)
{
	// mesh initialied in constructor
	// construct mesh from file if 2D, else construct based on n and p
	mesh->constructMesh(d == 2 ? fileNameNoExt : "");
	// allocate space for CSR stiffness matrix based on element DoFs
	mesh->allocateStiffness();
	// construct global stiffness matrix from local element stiffness matrices
	mesh->assembleStiffnessMatrix();
	// construct global load vector from local element load vectors
	mesh->assembleLoadVector(f);
	// apply boundary conditions,
	// set left and right boundaries in 1D,
	// or edge boundaries in 2D
	mesh->applyBoundaryConditions(u0, u1, boundary_u0, boundary_u1);

	// resize load and solution vectors based on number of nodes
	// since noNodes != n
	int nDoF = mesh->getNoNodes();
	globalLoad.resize(nDoF);
	solution.resize(nDoF);

	globalStiffness = mesh->stiffness;
	globalLoad = mesh->load;
	n = mesh->n;

	// initial guess for Picard iteration
	std::vector<double> U_n(nDoF, 0.0);
	std::vector<double> U_np1(nDoF, 0.0);

	for (int i=0; i<maxIter; i++)
	{
		// find a(U_n, v)
		std::vector<double> stiffnessProduct = dynamic_cast<FE_Mesh2D*>(mesh.get())->assembleStiffnessProduct(U_n);
		// find nonlinear term b(U_n; v)
		std::vector<double> load_NL = dynamic_cast<FE_Mesh2D*>(mesh.get())->assembleNonlinearLoad(U_n, q);

		// create rhs vector := (1 - alpha) * a(U_n, v) + alpha * [F - lambda * b(U_n; v)]
		// where F = globalLoad
		std::vector<double> rhs(nDoF, 0.0);
		for (int j=0; j<nDoF; j++)
		{
			rhs.at(j) = (1.0 - alpha) * stiffnessProduct.at(j) + alpha * (globalLoad.at(j) - lambda * load_NL.at(j));
		}

		// apply boundary conditions
		// TODO:

		// solve system A * u_n+1 = rhs
		Solver solver(n, p, globalStiffness, rhs, U_np1);
		solver.setupEigen();
		U_np1 = solver.solveEigen();

		// check for convergence through L^2 norm
		// ||u_n+1 - u_n||_{L^2} = sqrt((U_{n+1} - U_n)^T * globalStiffness * (U_{n+1} - U_n))
		std::vector<double> diff(nDoF, 0.0);
		for (int j=0; j<nDoF; j++)
		{
			diff[j] = U_np1[j] - U_n[j];
		}

		// find A * diff
		std::vector<double> Ax(nDoF, 0.0);
		for (int j=0; j<nDoF; j++)
		{
			for (int k=globalStiffness.row_start.at(j); k<globalStiffness.row_start.at(j+1); k++)
			{
				Ax[j] += globalStiffness.entries.at(k) * diff.at(globalStiffness.col_no.at(k));
			}
		}

		// find L2 norm
		double norm = 0.0;
		for (int j=0; j<nDoF; j++)
		{
			norm += Ax[j] * diff[j];
		}
		norm = pow(norm, 0.5);

		if (norm < tol)
		{
			std::cout << "Converged in " << i + 1 << " iterations." << std::endl;
			U_n = U_np1;
			break;
		}
		else if (i == maxIter - 1)
		{
			std::cout << "Did not converge in " << maxIter << " iterations." << std::endl;
			U_n = U_np1;
			break;
		}

		U_n = U_np1;
	}
	
	solution = U_n;
	return solution;
}

double FE_Solution::evaluateSolution(std::vector<double> x)
{
	// use mesh specific evaluation method
	return mesh->evaluateSolution(x, solution);
}

double FE_Solution::getMeshSize()
{
	// use mesh specific measurement
	return mesh->getMeshSize();
}

// dumps the triangulation and the solution on it, for post-processing
void FE_Solution::sendFieldToFile(int levels)
{
	FE_Mesh2D* mesh2D = dynamic_cast<FE_Mesh2D*>(mesh.get());

	if (mesh2D == nullptr)
	{
		std::cerr << "field output is two dimensional only" << std::endl;
		return;
	}

	mesh2D->sendFieldToFile(solution, levels);
}

void FE_Solution::evaluateDerivative(std::vector<double> x, double grad[2])
{
	// use mesh specific evaluation method
	mesh->evaluateDerivative(x, solution, grad);
}

double FE_Solution::getL2Error(double (*f_analytic)(const std::vector<double>&))
{
	double L2error = 0.0;

	if (d == 1)
	{
		GaussQuadrature quad;
		quad.assembleQuadrature(p+1);
		for (int i=0; i<n; i++)
		{
			// for each element get L and R nodes
			double nodeL = mesh->getNode(i)[0];
			double nodeR = mesh->getNode(i+1)[0];
			for (int j=0; j<p+1; j++)
			{
				double phi = 0.5 * (1.0 + quad.points[j]);
				double x = nodeL + phi * (nodeR - nodeL);
				double u_exact = f_analytic({x});
				double u_h = evaluateSolution({x});
	
				L2error += quad.weights[j] * 0.5 * (nodeR - nodeL) * pow(u_exact - u_h, 2);
			}
		}
	}
	else if (d == 2)
	{
		GaussQuadrature2D quad;
		quad.assembleQuadrature(p+1);
		for (int i=0; i<n; i++)
		{
			Element2D* elem = dynamic_cast<Element2D*>(mesh->elements[i].get());
			const std::vector<Point2D>& nodes = elem->nodes;

			double A[2][2];
			computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
			double detA = fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

			for (int j=0; j<static_cast<int>(quad.points.size()); j++)
			{
				Point2D point = mapToPhysical(nodes[0], nodes[1], nodes[2], quad.points[j]);
				double u_exact = f_analytic({point.x, point.y});
				double u_h = evaluateSolution({point.x, point.y});
				double diff = u_exact - u_h;
	
				L2error += quad.weights[j] * diff * diff * detA;
			}
		}
	}

	return pow(L2error, 0.5);
}

void FE_Solution::sendSolutionToFile(int noGridPoints, double (*f_analytic)(const std::vector<double>&))
{
	if (d == 1)
	{
		std::ofstream file("solution.csv");
		file << "x,u,u_exact\n";

		for (int i=0; i<noGridPoints + 1; i++)
		{
			double x = (double)i / (double)noGridPoints;
			file << x << ","
				 << evaluateSolution({x}) << ","
				 << f_analytic({x}) << "\n";
		}

		file.close();
	}
	else if (d == 2)
	{
		std::ofstream file("solution.csv");
		file << "x,y,u,u_exact\n";

		for (int i=0; i<noGridPoints + 1; i++)
		{
			for (int j=0; j<noGridPoints + 1; j++)
			{
				double x = (double)i / (double)noGridPoints;
				double y = (double)j / (double)noGridPoints;
				if (y < 0 && x > 0)
				{
					continue;
				}
				file << x << ","
					 << y << ","
					 << evaluateSolution({x, y}) << ","
					 << f_analytic({x, y}) << "\n";
			}
		}

		file.close();
	}
}

void FE_Solution::sendSolutionToFile(int noGridPoints)
{
	dynamic_cast<FE_Mesh2D*>(mesh.get())->sendSolutionToFile(noGridPoints, solution);
}