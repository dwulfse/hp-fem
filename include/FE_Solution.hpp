#ifndef FESOLUTIONHEADERDEF
#define FESOLUTIONHEADERDEF

#include <vector>
#include <fstream>
#include <memory>
#include "CSR_Matrix.hpp"
#include "FE_Mesh.hpp"

class FE_Solution
{
	public:
		int n; // number of elements
		int p; // polynomial degree
		int d; // dimension
		CSRMatrix globalStiffness;
		std::vector<double> globalLoad;
		std::vector<double> solution;
		// smart pointer uses automatic memory management
		// for added safety
		std::unique_ptr<FE_Mesh> mesh;

		// constructor
		FE_Solution(int n, int p, int d);

		// destructor
		~FE_Solution();
		
		std::vector<double> solve(double (*f)(const std::vector<double>&),
											  std::string fileNameNoExt="",
											  double u0=0.0, double u1=0.0,
											  bool boundary_u0=true, bool boundary_u1=true);
		std::vector<double> solveSemilinear(double (*f)(const std::vector<double>&),
											int q, double alpha, double lambda,
											int maxIter, double tol,
											std::string fileNameNoExt="",
											double u0=0.0, double u1=0.0,
											bool boundary_u0=true, bool boundary_u1=true);
		double evaluateSolution(std::vector<double> x);
		double getMeshSize();
		void evaluateDerivative(std::vector<double> x, double grad[2]);
		double getL2Error(double (*f_analytic)(const std::vector<double>&));
		void sendSolutionToFile(int n, double (*f_analytic)(const std::vector<double>&));
		void sendSolutionToFile(int noGridPoints);
};

#endif