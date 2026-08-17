#ifndef FEMESH2DHEADERDEF
#define FEMESH2DHEADERDEF

#include "FE_Mesh.hpp"
#include "Point2D.hpp"
#include <map>

class FE_Mesh2D : public FE_Mesh
{
	public:
		std::vector<Point2D> nodes;
		std::vector<bool> is_boundary;
		// edge DoFs
		// map holds (v1, v2): list of DoF indicies
		std::map<std::pair<int, int>, std::vector<int>> edge_dofs;

		// constructors
		FE_Mesh2D(int n, int p, int d);
		FE_Mesh2D();

		// destructor
		virtual ~FE_Mesh2D();

		// methods
		virtual std::vector<double> getNode(int i) override;
		virtual int getNoNodes() override;
		virtual double getMeshSize() override;
		virtual void constructMesh(std::string filename = "") override;
		virtual double evaluateSolution(std::vector<double> x, std::vector<double> solution) override;
		virtual void evaluateDerivative(std::vector<double> x, std::vector<double> solution, double grad[2]) override;
		virtual void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1) override;

		std::vector<double> assembleNonlinearLoad(const std::vector<double>& U, int q);
		std::vector<double> assembleStiffnessProduct(const std::vector<double>& U);
		void sendSolutionToFile(int n, const std::vector<double>& solution);
		void sendFieldToFile(const std::vector<double>& solution, int levels);
};

#endif