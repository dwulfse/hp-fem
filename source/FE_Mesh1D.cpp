#include "FE_Mesh1D.hpp"
#include "Element1D.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

// constructors
FE_Mesh1D::FE_Mesh1D()
{
}

FE_Mesh1D::FE_Mesh1D(int n, int p, int d)
 : FE_Mesh(n, p, d), nodes(p*n+1, 0.0)
{
}

// destructor
FE_Mesh1D::~FE_Mesh1D()
{
}

// get node i coordinates
std::vector<double> FE_Mesh1D::getNode(int i)
{
	return { nodes[i] };
}

// get number of nodes
int FE_Mesh1D::getNoNodes()
{
	return p*n + 1;
}

// construct mesh of n elements and n+1 nodes
// filename is unused: the 1D mesh is generated from n and p, not read from file
void FE_Mesh1D::constructMesh(std::string /*filename*/)
{
	const double h = 1.0 / n; // element size
	int int_i = n + 1; // start of extra nodes indexing

	// boundary nodes
	nodes.resize(n+1);
	nodes[0] = 0.0;
	for (int i=1; i<n+1; i++)
	{
		nodes[i] = i * h;
	}

	// elements
	std::vector<int> local_dofs(p+1);
	std::vector<double> element_nodes(2);
	for (int i=0; i<n; i++)
	{
		element_nodes[0] = nodes[i];
		element_nodes[1] = nodes[i+1];
		local_dofs[0] = i;
		local_dofs[1] = i+1;
		if (p > 1)
		{
			for (int j=2; j<p+1; j++)
			{
				local_dofs[j] = int_i++;
			}
		}
		elements[i] = std::make_unique<Element1D>(p, local_dofs, element_nodes);
	}
}

// evaluate solution at x
double FE_Mesh1D::evaluateSolution(std::vector<double> x, std::vector<double> solution)
{
	// find element containing x
	int elem = 0;
	for (int i=0; i<n; i++)
	{
		if (x[0] >= nodes[i] && x[0] <= nodes[i+1])
		{
			elem = i;
			break;
		}
	}

	// transform x to xi
	double xi = 2.0 * (x[0] - nodes[elem]) / (nodes[elem+1] - nodes[elem]) - 1.0;

	// evaluate u_h(x)
	double u_val = 0.0;
	const std::vector<int>& local_DoF = elements[elem]->local_DoF;
	for (int i=0; i<p+1; i++)
	{
		u_val += solution[local_DoF[i]] * elements[elem]->poly.evaluate_lobatto(i, xi);
	}

	return u_val;
}

// evaluate derivative at x
void FE_Mesh1D::evaluateDerivative(std::vector<double> x, std::vector<double> solution, double grad[2])
{
	// find element containing x
	int elem = 0;
	for (int i=0; i<n; i++)
	{
		if (x[0] >= nodes[i] && x[0] <= nodes[i+1])
		{
			elem = i;
			break;
		}
	}

	// transform x to xi
	double xi = 2.0 * (x[0] - nodes[elem]) / (nodes[elem+1] - nodes[elem]) - 1.0;

	// evaluate u_h(x)
	const std::vector<int>& local_DoF = elements[elem]->local_DoF;
	for (int i=0; i<p+1; i++)
	{
		grad[0] += solution[local_DoF[i]] * elements[elem]->poly.evaluate_lobatto_deriv(i, xi);
	}
}

// apply left and right boundary conditions
void FE_Mesh1D::applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	if (boundary_u0)
	{
		for (int i=stiffness.row_start[0]; i<stiffness.row_start[1]; i++)
		{
			stiffness.entries[i] = 0.0;
		}

		for (int i=1; i<static_cast<int>(stiffness.col_no.size()); i++)
		{
			if (stiffness.col_no[i] == 0)
			{
				stiffness.entries[i] = 0.0;
			}
		}

		stiffness(0, 0) = 1.0;
		load[0] = u0;
	}

	if (boundary_u1)
	{
		for (int i=stiffness.row_start[n]; i<stiffness.row_start[n+1]; i++)
		{
			stiffness.entries[i] = 0.0;
		}

		for (int i=0; i<static_cast<int>(stiffness.col_no.size()); i++)
		{
			if (stiffness.col_no[i] == n)
			{
				stiffness.entries[i] = 0.0;
			}
		}

		stiffness(n, n) = 1.0;
		load.at(n) = u1;
	}
}