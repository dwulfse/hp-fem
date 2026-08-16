#include "FE_Solution.hpp"
#include "Polynomial.hpp"
#include <iostream>
#include <string>
#include <cmath>

#include "Point2D.hpp"
#include "Element2D.hpp"
#include "PolynomialSpace.hpp"

int main()
{
	// setup problem
	const int polynomialDegree = 1; // polynomial degree
	const int dimension = 1; // dimension
	const bool semilinear = false; // semilinear problem

	std::vector<double> solution;

	// boundary conditions
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	if (dimension == 1)
	{
		const int n = 4; // number of elements
		// -u'' = f
		// f(x) = pi^2 sin(pi x)
		// u(x) = sin(pi x)
		auto f = [](const std::vector<double>& x) { return M_PI * M_PI * sin(M_PI * x[0]); }; // RHS function
		auto u_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]); }; // analytic solution
		// create and solve problem
		FE_Solution FEM(n, polynomialDegree, dimension);
		solution = FEM.solve(f, "", u0, u1, boundary_u0, boundary_u1);
		FEM.sendSolutionToFile(64, u_analytic);

		// output u vector
		std::cout << "u vector: " << std::endl;
		for (int i=0; i<static_cast<int>(solution.size()); i++)
		{
			std::cout << "u" << i << " = " << solution[i] << std::endl;
		}
		
		// output solution
		std::cout << "Solution:" << std::endl;
		for (int i=0; i<5; i++)
		{
			double x = i / 4.0;
			std::cout << "u_h(" << x << ") = " << FEM.evaluateSolution({x}) << std::endl;
		}
		// should be [0, 0.09375, 0.125, 0.09375, 0]

		// convergence analysis
		std::ofstream file("hp_error.csv");
		file << "p,h,error\n";

		for (int p=1; p<=3; p++)
		{
			for (int n=10; n<=20; n++)
			{
				FE_Solution FEM(n, p, dimension);
				std::vector<double> solution = FEM.solve(f, "", u0, u1, boundary_u0, boundary_u1);
				double error = FEM.getL2Error(u_analytic);
				double h = FEM.getMeshSize();
				file << p << "," << h << "," << error << "\n";
			}
		}

		file.close();
	}
	else if (dimension == 2)
	{
		if (!(semilinear))
		{
			// -lap u = f
			// f(x, y) = 2pi^2 sin(pi x) sin(pi y)
			// u(x, y) = sin(pi x) sin(pi y)
			auto f = [](const std::vector<double>& x) { return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // RHS function
			auto u_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // analytic solution
			std::string filenameNoExt = "6";
			// create and solve problem
			FE_Solution FEM(1, polynomialDegree, dimension);
			solution = FEM.solve(f, filenameNoExt, u0, u1, boundary_u0, boundary_u1);
			FEM.sendSolutionToFile(32, u_analytic);

			// output u vector
			std::cout << "u vector: " << std::endl;
			for (int i=0; i<static_cast<int>(solution.size()); i++)
			{
				std::cout << "u" << i << " = " << solution[i] << std::endl;
			}

			// convergence analysis
			// domain.9 is excluded: it holds 15565 elements against domain.8's 512,
			// so it is not a uniform refinement of the same family and would skew the study
			std::vector<std::string> mesh_files = {"4", "5", "6", "7", "8"};
			std::ofstream file("hp_error.csv");
			file << "p,h,error\n";

			for (int p=1; p<=3; p++)
			{
				for (int i=0; i<static_cast<int>(mesh_files.size()); i++)
				{
					FE_Solution FEM(1, p, dimension);
					FEM.solve(f, mesh_files[i]);
					double error = FEM.getL2Error(u_analytic);
					// h is measured from the mesh, not inferred from the filename
					double h = FEM.getMeshSize();
					file << p << "," << h << "," << error << "\n";
					std::cout << "p = " << p << ", h = " << h << ", error = " << error << std::endl;
				}
			}

			file.close();
		}
		else
		{
			// -lap u + lambda u^{2q+1}= f
			// f(x, y) = 2pi^2 sin(pi x) sin(pi y)
			auto f = [](const std::vector<double>& x) { return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // RHS function
			std::string filenameNoExt = "6";

			int q = 100;
			double alpha = 0.7;
			double lambda = 100.0;
			int maxIter = 100;
			double tol = 1e-8;

			// create and solve problem
			FE_Solution FEM(1, polynomialDegree, 2);

			solution = FEM.solveSemilinear(f, q, alpha, lambda, maxIter, tol, filenameNoExt);
			FEM.sendSolutionToFile(32, f);

			// output u vector
			std::cout << "u vector: " << std::endl;
			for (int i=0; i<static_cast<int>(solution.size()); i++)
			{
				std::cout << "u" << i << " = " << solution[i] << std::endl;
			}
		}
	}

	return 0;
}