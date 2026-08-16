#ifndef FEMESH1DHEADERDEF
#define FEMESH1DHEADERDEF

#include "FE_Mesh.hpp"

class FE_Mesh1D : public FE_Mesh
{
	public:
		std::vector<double> nodes;

		// constructors
		FE_Mesh1D(int n, int p, int d);
		FE_Mesh1D();

		// destructor
		virtual ~FE_Mesh1D();

		// methods
		virtual std::vector<double> getNode(int i) override;
		virtual int getNoNodes() override;
		virtual double getMeshSize() override;
		virtual void constructMesh(std::string filename = "") override;
		virtual double evaluateSolution(std::vector<double> x, std::vector<double> solution) override;
		virtual void evaluateDerivative(std::vector<double> x, std::vector<double> solution, double grad[2]) override;
		virtual void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1) override;
};

#endif