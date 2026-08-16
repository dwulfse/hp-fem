#include "Element2D.hpp"
#include "GaussQuadrature2D.hpp"
#include "Helper.hpp"
#include <cmath>
#include <iostream>

// constructors
Element2D::Element2D()
{
}

Element2D::Element2D(int id, int p, std::vector<int>& local_DoF, std::vector<Point2D>& nodes)
 : Element(id, p, local_DoF), nodes(nodes)
{
}

// destructor
Element2D::~Element2D()
{
}

std::vector<std::vector<double>> Element2D::getLocalStiffness()
{
	// find degrees of freedom and allocate local stiffness matrix
	int nDoF = (p+1) * (p+2) / 2;
	std::vector<std::vector<double>> stiffness(nDoF, std::vector<double>(nDoF, 0.0));

	// create quadrature
	GaussQuadrature2D quad;
	quad.assembleQuadrature(p+1);

	// find affine transformation matrix and inverse
	double A[2][2];
	computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
	double detA = fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

	double A_inv[2][2];
	invertAffineMatrix(A, A_inv);

	// compute local stiffness matrix, integrated by quadrature
	for (int k=0; k<static_cast<int>(quad.points.size()); k++)
	{
		std::vector<std::vector<double>> gradRef(nDoF, std::vector<double>(2, 0.0));
		for (int i=0; i<nDoF; i++)
		{
			double gradRef_i[2] = {0.0, 0.0};
			// evaluate ith basis function at quadrature point
			poly.basis_2D_grad(i, quad.points[k].x, quad.points[k].y, gradRef_i);
			gradRef[i][0] = gradRef_i[0];
			gradRef[i][1] = gradRef_i[1];
		}
		
		// transform reference gradients to physical gradients
		// grad_phys phi = A^{-T} grad_ref phi
		std::vector<std::vector<double>> gradPhys(nDoF, std::vector<double>(2, 0.0));
		for (int i=0; i<nDoF; i++)
		{
			gradPhys[i][0] = A_inv[0][0] * gradRef[i][0] + A_inv[1][0] * gradRef[i][1];
			gradPhys[i][1] = A_inv[0][1] * gradRef[i][0] + A_inv[1][1] * gradRef[i][1];
		}

		for (int i=0; i<nDoF; i++)
		{
			for (int j=0; j<nDoF; j++)
			{
				// scaling factor is |T_phys| / |T_ref| = (detA / 2) / 2 = detJ / 4
				stiffness[i][j] += (detA) * quad.weights[k] *
						(gradPhys[i][0] * gradPhys[j][0] + gradPhys[i][1] * gradPhys[j][1]);
			}
		}
	}

	return stiffness;
}

std::vector<double> Element2D::getLocalStiffnessProduct(const std::vector<double>& U)
{
	// find degrees of freedom and allocate local stiffness matrix
	int nDoF = (p+1) * (p+2) / 2;
	std::vector<double> stiffnessProduct(nDoF, 0.0);

	// create quadrature
	GaussQuadrature2D quad;
	quad.assembleQuadrature(p+1);

	// find affine transformation matrix and inverse
	double A[2][2];
	computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
	double detA = fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

	double A_inv[2][2];
	invertAffineMatrix(A, A_inv);

	// compute local stiffness matrix, integrated by quadrature
	for (int k=0; k<static_cast<int>(quad.points.size()); k++)
	{
		std::vector<std::vector<double>> gradRef(nDoF, std::vector<double>(2, 0.0));
		for (int i=0; i<nDoF; i++)
		{
			double gradRef_i[2] = {0.0, 0.0};
			poly.basis_2D_grad(i, quad.points[k].x, quad.points[k].y, gradRef_i);
			gradRef.at(i)[0] = gradRef_i[0];
			gradRef.at(i)[1] = gradRef_i[1];
		}

		std::vector<std::vector<double>> gradPhys(nDoF, std::vector<double>(2, 0.0));
		for (int i=0; i<nDoF; i++)
		{
			gradPhys[i][0] = A_inv[0][0] * gradRef[i][0] + A_inv[1][0] * gradRef[i][1];
			gradPhys[i][1] = A_inv[0][1] * gradRef[i][0] + A_inv[1][1] * gradRef[i][1];
		}

		double grad_U[2] = {0.0, 0.0};
		for (int j=0; j<nDoF; j++)
		{
			grad_U[0] += U.at(local_DoF[j]) * gradPhys[j][0];
			grad_U[1] += U.at(local_DoF[j]) * gradPhys[j][1];
		}

		for (int i=0; i<nDoF; i++)
		{
			// grad U dot grad phi
			double integrand = (grad_U[0] * gradPhys[i][0] + grad_U[1] * gradPhys[i][1]);
			stiffnessProduct[i] += (detA) * quad.weights[k] * integrand;
		}
	}

	return stiffnessProduct;
}

std::vector<double> Element2D::getLocalLoad(double (*f)(const std::vector<double>&))
{
	// find degrees of freedom and allocate local load vector
	int nDoF = local_DoF.size();		// (p+1) * (p+2) / 2;
	std::vector<double> load(nDoF, 0.0);

	// create quadrature
	GaussQuadrature2D quad;
	quad.assembleQuadrature(p+1);

	// find affine transformation matrix and inverse
	double A[2][2];
	computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
	double detA = fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

	// compute local load vector, integrated by quadrature
	for (int k=0; k<static_cast<int>(quad.points.size()); k++)
	{
		Point2D point = mapToPhysical(nodes[0], nodes[1], nodes[2], quad.points[k]);
		double fx = f({point.x, point.y});

		for (int i=0; i<nDoF; i++)
		{
			double phi = poly.basis_2D(i, quad.points[k].x, quad.points[k].y);
			load[i] += (detA / 4.0) * quad.weights[k] * fx * phi;
		}
	}

	return load;
}

std::vector<double> Element2D::getLocalNonlinearLoad(const std::vector<double>& U, int q)
{
	// find degrees of freedom and allocate local load vector
	int nDoF = local_DoF.size();
	std::vector<double> loadNL(nDoF, 0.0);

	// create quadrature
	GaussQuadrature2D quad;
	quad.assembleQuadrature(p+1);

	// find affine transformation matrix and inverse
	double A[2][2];
	computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
	double detA = fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

	for (int k=0; k<static_cast<int>(quad.points.size()); k++)
	{
		double xi1 = quad.points[k].x;
		double xi2 = quad.points[k].y;

		// reconstruct U at quadrature point
		double u_val = 0.0;
		for (int j=0; j<nDoF; j++)
		{
			double phi = poly.basis_2D(j, xi1, xi2);
			u_val += U.at(local_DoF[j]) * phi;
		}

		//compute nonlinear U^{2q+1} value
		double u_pow = pow(u_val, 2*q+1);

		// assemble contributions for each local basis function
		for (int i=0; i<nDoF; i++)
		{
			double phi = poly.basis_2D(i, xi1, xi2);
			loadNL[i] += (detA / 4.0) * quad.weights[k] * u_pow * phi;
		}
	}
	return loadNL;
}