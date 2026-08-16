#include "Element1D.hpp"
#include "GaussQuadrature.hpp"

// constructors
Element1D::Element1D()
{
}

Element1D::Element1D(int p, std::vector<int>& local_DoF, std::vector<double>& nodes)
 : Element(p, local_DoF), nodes(nodes)
{
}

// destructor
Element1D::~Element1D()
{
}

std::vector<std::vector<double>> Element1D::getReferenceStiffness()
{
	std::vector<std::vector<double>> refStiffness(p+1, std::vector<double>(p+1, 0.0));
	GaussQuadrature quad;
	quad.assembleQuadrature(p+1);
	
	double gauss;
	for (int k=0; k<p+1; k++)
	{
		gauss = quad.points[k];
		for (int i=0; i<p+1; i++)
		{
			for (int j=0; j<p+1; j++)
			{
				refStiffness[i][j] += quad.weights[k] * poly.evaluate_lobatto_deriv(i, gauss) * poly.evaluate_lobatto_deriv(j, gauss);
			}
		}
	}

	return refStiffness;
}

std::vector<std::vector<double>> Element1D::getLocalStiffness()
{
	std::vector<std::vector<double>> stiffness(p+1, std::vector<double>(p+1, 0.0));
	std::vector<std::vector<double>> refStiffness = getReferenceStiffness();
	double J = (nodes.back() - nodes.front()) / 2.0;

	for (int i=0; i<p+1; i++)
	{
		for (int j=0; j<p+1; j++)
		{
			stiffness[i][j] = refStiffness[i][j] / J;
		}
	}

	return stiffness;
}

std::vector<double> Element1D::getLocalLoad(double (*f)(const std::vector<double>&))
{
	std::vector<double> local_load(p+1, 0.0);
	GaussQuadrature quad;
	quad.assembleQuadrature(p+1);

	double h = nodes.back() - nodes.front();
	double J = h / 2.0;

	double gauss;
	for (int k=0; k<static_cast<int>(quad.points.size()); k++)
	{
		gauss = quad.points[k];
		double x = 0.5 * h * gauss + 0.5 * (nodes.back() + nodes.front());
		for (int i=0; i<p+1; i++)
		{
			local_load[i] += quad.weights[k] * J * f({x}) * poly.evaluate_lobatto(i, gauss);
		}
	}

	return local_load;
}