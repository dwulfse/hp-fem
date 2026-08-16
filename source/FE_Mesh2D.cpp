#include "FE_Mesh2D.hpp"
#include "Element2D.hpp"
#include "GaussQuadrature2D.hpp"
#include "Helper.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

// make pairs of ordered edges based on two vertices
std::pair<int, int> makeEdgePair(int i, int j)
{
	return i < j ? std::make_pair(i, j) : std::make_pair(j, i);
}

// constructors
FE_Mesh2D::FE_Mesh2D()
{
}

FE_Mesh2D::FE_Mesh2D(int n, int p, int d)
 : FE_Mesh(n, p, d), nodes(0), is_boundary(0), edge_dofs()
{
}

// destructor
FE_Mesh2D::~FE_Mesh2D()
{
}

// get node i coordinates
std::vector<double> FE_Mesh2D::getNode(int i)
{
	return { nodes[i].x, nodes[i].y };
}

// get number of DoFs
int FE_Mesh2D::getNoNodes()
{
	// no. DoFs = [no. vertices] + (p - 1) [no. edges] + (p - 1)(p - 2)/2 [no. bubbles]
	return nodes.size() +
		(p - 1) * edge_dofs.size() +
		(p - 1) * (p - 2) / 2 * n;
}

// load mesh from Triangle generated files
void FE_Mesh2D::constructMesh(std::string filename)
{
	// open files
	std::ifstream eleFile(std::string("domain.") + filename + ".ele");
	std::ifstream nodeFile(std::string("domain.") + filename + ".node");

	if (!eleFile.good() || !nodeFile.good())
	{
		throw std::runtime_error("Mesh files not found");
	}

	std::string line;

	// read node header
	std::getline(nodeFile, line);
	if (line.empty())
	{
		throw std::runtime_error("Empty node file header");
	}
	std::istringstream nodeHeader(line);
	int noVertices, dimension, noAttributesNode, boundary;
	nodeHeader >> noVertices >> dimension >> noAttributesNode >> boundary;

	// resize vectors to number of vertices
	nodes.resize(noVertices);
	load.resize(noVertices, 0.0);
	is_boundary.resize(noVertices);

	// read nodes
	while (std::getline(nodeFile, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream iss(line);
		Point2D node;
		int boundary_flag;
		// line format: id x y [boundary_flag]
		if (!(iss >> node.id >> node.x >> node.y >> boundary_flag))
		{
			std::cerr << "Couldn't read node file line: " << line << std::endl;
			continue;
		}
		node.id--; // 0-indexed
		if (node.id < 0 || node.id >= noVertices)
		{
			std::cerr << "Node id " << node.id << " out of range" << std::endl;
			continue;
		}
		nodes[node.id] = node;
		if (boundary_flag == 1)
		{
			is_boundary[node.id] = true;
		}
		else
		{
			is_boundary[node.id] = false;
		}
	}

	// read element header
	std::getline(eleFile, line);
	if (line.empty())
	{
		throw std::runtime_error("Empty element file header");
	}
	std::istringstream eleHeader(line);
	int noElems, noNodesElem, noAttributes;
	eleHeader >> noElems >> noNodesElem >> noAttributes;
	if (noNodesElem != 3)
	{
		throw std::runtime_error("Only triangles supported");
	}

	// in 2D, n (number of elements) read from file
	// resize elements vector
	n = noElems;
	elements.resize(noElems);
	is_boundary.resize(noVertices + (p-1) * (3 * n + noVertices - 3) / 2 + (p-1) * (p-2) * n / 2, false);

	int DoFi = noVertices; // start of extra DoFs indexing

	int nDoFs = (p+1) * (p+2) / 2; // number of DoFs per element
	int nBubbleDoFs = (p-1) * (p-2) / 2; // number of bubble DoFs per element

	// read elements
	// setup element nodes and local DoFs to populate for each element
	// local DoFs will hold vertex (connecivity), edge, and bubble DoFs
	std::vector<int> local_dofs(nDoFs);
	std::vector<Point2D> element_nodes(noNodesElem);
	int id;
	while(std::getline(eleFile, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream iss(line);
		// line format: id n1 n2 n3 [attribute]
		// read in element id
		if (!(iss >> id))
		{
			std::cerr << "Couldn't read element id from line" << line << std::endl;
			continue;
		}
		id--; // 0-indexed
		if (id < 0 || id >= noElems)
		{
			std::cerr << "Element id " << id << " out of range" << std::endl;
			continue;
		}
		// noNodesElem = 3 for triangular elements
		// read connectivity
		double connectivity[3];
		for (int i=0; i<noNodesElem; i++)
		{
			if (!(iss >> connectivity[i]))
			{
				std::cerr << "Couldn't read local dof for element " << id+1 << " from line" << line << std::endl;
				break;
			}
			connectivity[i]--; // 0-indexed
			if (connectivity[i] < 0 || connectivity[i] >= noVertices)
			{
				std::cerr << "Local dof " << connectivity[i] << " out of range" << std::endl;
				break;
			}
			element_nodes[i] = nodes[connectivity[i]];
		}

		// add vertex DoFs
		for (int i=0; i<noNodesElem; i++)
		{
			local_dofs[i] = connectivity[i];
		}

		// add edge DoFs
		for (int i=0; i<noNodesElem; i++)
		{
			int j = (i+1) % noNodesElem;
			int v1 = local_dofs[i], v2 = local_dofs[j];

			std::pair<int, int> edge = makeEdgePair(v1, v2);
			std::vector<int> dof(p-1);
			if (edge_dofs.find(edge) == edge_dofs.end())
			{
				dof.resize(p - 1);
				for (int k=0; k<p-1; k++)
				{
					is_boundary[DoFi] = is_boundary.at(v1) && is_boundary.at(v2);
					dof[k] = DoFi++;
				}
				edge_dofs[edge] = dof;
			}
			else
			{
				dof = edge_dofs[edge];
			}

			int edge_i = 3 + i * (p - 1);
			for (int k=0; k<p-1; k++)
			{
				local_dofs[edge_i + k] = dof[k];
			}
		}

		// add bubble DoFs
		int bubble_i = 3 + 3 * (p - 1);
		for (int i=0; i<nBubbleDoFs; i++)
		{
			is_boundary[DoFi] = false;
			local_dofs[bubble_i + i] = DoFi++;
		}
		// set Element to be type Element2D with read in DoFs and nodes
		elements[id] = std::make_unique<Element2D>(id, p, local_dofs, element_nodes);
	}

	// close files
	eleFile.close();
	nodeFile.close();
}

// evaluate solution at x
double FE_Mesh2D::evaluateSolution(std::vector<double> x, std::vector<double> solution)
{
	for (int i=0; i<static_cast<int>(elements.size()); i++)
	{
		Element2D* elem = dynamic_cast<Element2D*>(elements[i].get());
		PolynomialSpace poly = elem->poly;
		Point2D P = elem->nodes[0], Q = elem->nodes[1], R = elem->nodes[2];

		double A[2][2];
		computeAffineMatrix(P, Q, R, A);

		double A_inv[2][2];
		if (!invertAffineMatrix(A, A_inv))
		{
			std::runtime_error("Singular transformation matrix");
			continue;
		}

		double xi1 = 2.0 * (A_inv[0][0] * (x[0] - P.x) + A_inv[0][1] * (x[1] - P.y)) - 1.0;
		double xi2 = 2.0 * (A_inv[1][0] * (x[0] - P.x) + A_inv[1][1] * (x[1] - P.y)) - 1.0;

		if (xi1 < -1.0 || xi1 > 1.0 || xi2 < -1.0 || xi2 > 1.0)
		{
			continue;
		}

		double lambda[3];
		poly.evaluate_affine(xi1, xi2, lambda);

		if (lambda[0] >= -1e-12 && lambda[1] >= -1e-12 && lambda[2] >= -1e-12)
		{
			const std::vector<int>& elem_dof = elem->local_DoF;

			double u_val = 0.0;
			for (int j=0; j<static_cast<int>(elem_dof.size()); j++)
			{
				double phi = poly.basis_2D(j, xi1, xi2);
				u_val += solution[elem_dof[j]] * phi;
			}
			return u_val;
		}
	}
	return 0.0;
}

// evaluate derivative at x
void FE_Mesh2D::evaluateDerivative(std::vector<double> x, std::vector<double> solution, double grad[2])
{
	for (int i=0; i<static_cast<int>(elements.size()); i++)
	{
		Element2D* elem = dynamic_cast<Element2D*>(elements[i].get());
		PolynomialSpace poly = elem->poly;
		Point2D P = elem->nodes[0], Q = elem->nodes[1], R = elem->nodes[2];

		double A[2][2];
		computeAffineMatrix(P, Q, R, A);

		double A_inv[2][2];
		if (!invertAffineMatrix(A, A_inv))
		{
			std::runtime_error("Singular transformation matrix");
			continue;
		}

		double xi1 = 2.0 * (A_inv[0][0] * (x[0] - P.x) + A_inv[0][1] * (x[1] - P.y)) - 1.0;
		double xi2 = 2.0 * (A_inv[1][0] * (x[0] - P.x) + A_inv[1][1] * (x[1] - P.y)) - 1.0;

		if (xi1 < -1.0 || xi1 > 1.0 || xi2 < -1.0 || xi2 > 1.0)
		{
			continue;
		}

		double lambda[3];
		poly.evaluate_affine(xi1, xi2, lambda);

		if (lambda[0] >= -1e-12 && lambda[1] >= -1e-12 && lambda[2] >= -1e-12)
		{
			const std::vector<int>& elem_dof = elem->local_DoF;

			for (int i=0; i<static_cast<int>(elem_dof.size()); i++)
			{
				double phi_grad[2];

				poly.basis_2D_grad(i, xi1, xi2, phi_grad);
				grad[0] += solution[elem_dof[i]] * phi_grad[0];
				grad[1] += solution[elem_dof[i]] * phi_grad[1];
			}
		}
	}
}

// send solution to file
void FE_Mesh2D::sendSolutionToFile(int noGridPoints, const std::vector<double>& solution)
{
	std::ofstream file("solution.csv");
	file << "x,y,u\n";

	GaussQuadrature2D quad;
	quad.assembleQuadrature(noGridPoints);

	for (int i=0; i<static_cast<int>(elements.size()); i++)
	{
		Element2D* elem = dynamic_cast<Element2D*>(elements[i].get());
		const std::vector<Point2D>& nodes = elem->nodes;
		for (int j=0; j<noGridPoints; j++)
		{
			Point2D point = mapToPhysical(nodes[0], nodes[1], nodes[2], quad.points[j]);

			double u_val = 0.0;
			for (int k=0; k<static_cast<int>(elem->local_DoF.size()); k++)
			{
				double phi = elem->poly.basis_2D(k, quad.points[j].x, quad.points[j].y);
				u_val += solution[elem->local_DoF[k]] * phi;
			}

			file << point.x << "," << point.y << "," << u_val << "\n";
		}
	}

	file.close();
}

// apply boundary conditions
void FE_Mesh2D::applyBoundaryConditions(double u_val, double /*unused*/, bool apply_boundary, bool /*unused*/)
{
	if (!apply_boundary) return;

	for (int k=0; k<static_cast<int>(is_boundary.size()); k++)
	{
		if (is_boundary[k])
		{
			// zero out row
			for (int i=stiffness.row_start[k]; i<stiffness.row_start[k+1]; i++)
			{
				stiffness.entries[i] = 0.0;
			}

			// zero out column
			for (int i=0; i<static_cast<int>(stiffness.row_start.size()); i++)
			{
				for (int j=stiffness.row_start[i]; j<stiffness.row_start[i+1]; j++)
				{
					if (stiffness.col_no[j] == k)
					{
						stiffness.entries[j] = 0.0;
					}
				}
			}

			// set diagonal to 1
			stiffness(k, k) = 1.0;
			// set load vector to boundary value
			load[k] = u_val;
		}
	}
}

// assemble nonlinear load vector
std::vector<double> FE_Mesh2D::assembleNonlinearLoad(const std::vector<double>& U, int q)
{
	int totalDoFs = getNoNodes();
	std::vector<double> loadNL(totalDoFs, 0.0);

	for (int i=0; i<static_cast<int>(elements.size()); i++)
	{
		Element2D* elem = dynamic_cast<Element2D*>(elements[i].get());
		std::vector<double> local_loadNL = elem->getLocalNonlinearLoad(U, q);
		const std::vector<int>& elem_dof = elem->local_DoF;

		for (int j=0; j<static_cast<int>(elem_dof.size()); j++)
		{
			loadNL.at(elem_dof.at(j)) += local_loadNL.at(j);
		}
	}

	// apply boundary conditions
	for (int i=0; i<totalDoFs; i++)
	{
		if (is_boundary[i])
		{
			loadNL[i] = 0.0;
		}
	}

	return loadNL;
}

// assemble stiffness matrix for given U
std::vector<double> FE_Mesh2D::assembleStiffnessProduct(const std::vector<double>& U)
{
	int totalDoFs = getNoNodes();
	std::vector<double> stiffnessProduct(totalDoFs, 0.0);

	for (int i=0; i<static_cast<int>(elements.size()); i++)
	{
		Element2D* elem = dynamic_cast<Element2D*>(elements[i].get());
		std::vector<double> local_stiffness = elem->getLocalStiffnessProduct(U);
		const std::vector<int>& elem_dof = elem->local_DoF;

		for (int j=0; j<static_cast<int>(elem_dof.size()); j++)
		{
			stiffnessProduct.at(elem_dof.at(j)) += local_stiffness.at(j);
		}
	}

	// apply boundary conditions
	for (int i=0; i<totalDoFs; i++)
	{
		if (is_boundary[i])
		{
			stiffnessProduct[i] = 0.0;
		}
	}

	return stiffnessProduct;
}
