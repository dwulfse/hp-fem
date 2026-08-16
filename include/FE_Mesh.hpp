#ifndef FEMESHHEADERDEF
#define FEMESHHEADERDEF

#include <vector>
#include <string>
#include <memory>
#include "Element.hpp"
#include "CSR_Matrix.hpp"

class FE_Mesh
{
	public:
		int n; // number of elements
		int p; // polynomial degree
		int d; // dimension
		std::vector<std::unique_ptr<Element>> elements;
		CSRMatrix stiffness;
		std::vector<double> load;

		// constructors
		// default constructor
		FE_Mesh() {};
		// constructor with n and p for 1D
		FE_Mesh(int n, int p, int d) : n(n), p(p), d(d), elements(n), stiffness(), load(p*n+1, 0.0) {};

		// destructor
		virtual ~FE_Mesh() {};

		// pure virtual method
		virtual std::vector<double> getNode(int i) = 0;
		virtual int getNoNodes() = 0;
		// mesh size h, the largest element diameter, for convergence studies
		virtual double getMeshSize() = 0;
		virtual void constructMesh(std::string filename = "") = 0;
		virtual double evaluateSolution(std::vector<double> x, std::vector<double> solution) = 0;
		virtual void evaluateDerivative(std::vector<double> x, std::vector<double> solution, double grad[2]) = 0;
		virtual void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1) = 0;

		// virtual methods
		virtual void allocateStiffness();
		virtual void assembleStiffnessMatrix();
		virtual void assembleLoadVector(double (*f)(const std::vector<double>&));
};

#endif