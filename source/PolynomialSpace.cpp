#include "PolynomialSpace.hpp"
#include "Polynomial.hpp"

#include <cmath>
#include <stdexcept>

// constructors
PolynomialSpace::PolynomialSpace()
{
}

PolynomialSpace::PolynomialSpace(int p)
 : p(p)
{
}

// destructor
PolynomialSpace::~PolynomialSpace()
{
}

double PolynomialSpace::evaluate(int i, double x)
{
	if (i == 0)
	{
		return (1.0 - x) / 2.0;
	}
	else if (i == 1)
	{
		return (1.0 + x) / 2.0;
	}
	return 0.0;
}

// x is unused: the vertex shape functions are linear, so their derivatives are constant
double PolynomialSpace::evaluate_deriv(int i, double /*x*/)
{
	if (i == 0)
	{
		return 1.0 / 2.0;
	}
	else if (i == 1)
	{
		return -1.0 / 2.0;
	}
	return 0.0;
}

double PolynomialSpace::evaluate_lobatto(int i, double x)
{
	return lobatto(i, x);
}

double PolynomialSpace::evaluate_lobatto_deriv(int i, double x)
{
	if (i == 0)
	{
		return -0.5;
	}
	else if (i == 1)
	{
		return 0.5;
	}
	else
	{
		return (legendre_deriv(i-2, x) - legendre_deriv(i, x)) / sqrt(2.0 * (2.0 * i - 1.0));
	}
}

void PolynomialSpace::evaluate_affine(double xi1, double xi2, double lambda[3])
{
	lambda[0] = (xi2 + 1.0) / 2.0;
	lambda[1] = - (xi1 + xi2) / 2.0;
	lambda[2] = (xi1 + 1.0) / 2.0;
}

// xi1, xi2 are unused: the affine coordinates are linear, so their gradients are constant
void PolynomialSpace::evaluate_affine_grad(double /*xi1*/, double /*xi2*/, double grad[3][2])
{
	grad[0][0] = 0.0;
	grad[0][1] = 0.5;
	grad[1][0] = -0.5;
	grad[1][1] = -0.5;
	grad[2][0] = 0.5;
	grad[2][1] = 0.0;
}

double PolynomialSpace::kernel(int k_2, double x)
{
	// phi_{k-2}(x) is st. l_k(x) = l_0(x) * l_1(x) * phi_{k-2}(x)
	double l0 = (1.0 - x) / 2.0;
	double l1 = (1.0 + x) / 2.0;
	double lk = lobatto(k_2 + 2, x);
	return lk / (l0 * l1);
}

double PolynomialSpace::kernel_deriv(int k_2, double x)
{
	double l0 = (1.0 - x) / 2.0;
	double l1 = (1.0 + x) / 2.0;
	double lk = lobatto(k_2 + 2, x);
	double dlk = evaluate_lobatto_deriv(k_2 + 2, x);
	double denom = l0 * l1;
	double dDenom = - x / 2.0;
	return (dlk * denom - lk * dDenom) / (denom * denom);
}

double PolynomialSpace::evaluate_edge(int k, double lambda_a, double lambda_b)
{
	if (fabs(lambda_a) < 1e-12 || fabs(lambda_b) < 1e-12)
	{
		return 0.0;
	}
	return lambda_a * lambda_b * kernel(k-2, lambda_b - lambda_a);
}

double PolynomialSpace::evaluate_bubble(int n1, int n2, const double lambda[3])
{
	if (fabs(lambda[0]) < 1e-12 || fabs(lambda[1]) < 1e-12 || fabs(lambda[2]) < 1e-12)
	{
		return 0.0;
	}
	return lambda[0] * lambda[1] * lambda[2] * kernel(n1-1, lambda[2] - lambda[1]) * kernel(n2-1, lambda[1] - lambda[0]);
}

double PolynomialSpace::basis_2D(int i, double x1, double x2)
{
	// validate input
	if (i > (p+1)*(p+2) / 2)
	{
		throw std::runtime_error("Invalid basis function index");
	}

	// get affine coordinates
	double lambda[3];
	evaluate_affine(x1, x2, lambda);

	// vertex functions
	if (i == 0) return lambda[1];
	if (i == 1) return lambda[2];
	if (i == 2) return lambda[0];

	// edge functions
	int noVertices = 3;
	int modesPEdge = p - 1;
	int noEdgeFuncs = 3 * modesPEdge;
	if (i < noVertices + noEdgeFuncs)
	{
		int edge_i = (i - noVertices) / modesPEdge;
		// edge modes run k = 2, ..., p; evaluate_edge takes k and uses kernel k-2
		int k = ((i - noVertices) % modesPEdge) + 2;

		switch(edge_i)
		{
			case 0:
				return evaluate_edge(k, lambda[1], lambda[2]);
			case 1:
				return evaluate_edge(k, lambda[2], lambda[0]);
			case 2:
				return evaluate_edge(k, lambda[0], lambda[1]);
			default:
				throw std::runtime_error("Invalid edge index");
		}
	}

	// bubble functions
	int bubble_i = i - noVertices - noEdgeFuncs;

	// find n_1, n_2 corresponding to index
	int n1 = 1, n2 = 1, counter = 0;
	bool found = false;
	for (n1=1; n1<p-1; n1++)
	{
		for (n2=1; n2<p-n1; n2++)
		{
			if (counter == bubble_i)
			{
				found = true;
				break;
			}
			counter++;
		}
		if (found) break;
	}

	return evaluate_bubble(n1, n2, lambda);
}

void PolynomialSpace::basis_2D_grad(int i, double x1, double x2, double grad[2])
{
	// validate input
	if (i > (p+1)*(p+2) / 2)
	{
		throw std::runtime_error("Invalid basis function index");
	}

	// get affine coordinates
	double lambda[3];
	evaluate_affine(x1, x2, lambda);

	// get affine derivatives
	double grad_lambda[3][2];
	evaluate_affine_grad(x1, x2, grad_lambda);

	// vertex functions
	if (i == 0)
	{
		grad[0] = grad_lambda[1][0];
		grad[1] = grad_lambda[1][1];
		return;
	}
	if (i == 1)
	{
		grad[0] = grad_lambda[2][0];
		grad[1] = grad_lambda[2][1];
		return;
	}
	if (i == 2)
	{
		grad[0] = grad_lambda[0][0];
		grad[1] = grad_lambda[0][1];
		return;
	}

	// edge functions
	int noVertices = 3;
	int modesPEdge = p - 1;
	int noEdgeFuncs = 3 * modesPEdge;
	if (i < noVertices + noEdgeFuncs)
	{
		int edge_i = (i - noVertices) / modesPEdge;
		// edge modes run k = 2, ..., p; the kernel below is indexed k-2
		int k = ((i - noVertices) % modesPEdge) + 2;

		double product, diff;
		double dProd[2], dDiff[2];

		switch (edge_i)
		{
			case 0:
				product = lambda[1] * lambda[2];
				diff = lambda[2] - lambda[1];
				// (lambda_1 * lambda_2)' = lambda_1' * lambda_2 + lambda_1 * lambda_2'
				dProd[0] = lambda[2] * grad_lambda[1][0] + lambda[1] * grad_lambda[2][0];
				dProd[1] = lambda[2] * grad_lambda[1][1] + lambda[1] * grad_lambda[2][1];
				dDiff[0] = grad_lambda[2][0] - grad_lambda[1][0];
				dDiff[1] = grad_lambda[2][1] - grad_lambda[1][1];
				break;
			case 1:
				product = lambda[2] * lambda[0];
				diff = lambda[0] - lambda[2];
				dProd[0] = lambda[0] * grad_lambda[2][0] + lambda[2] * grad_lambda[0][0];
				dProd[1] = lambda[0] * grad_lambda[2][1] + lambda[2] * grad_lambda[0][1];
				dDiff[0] = grad_lambda[0][0] - grad_lambda[2][0];
				dDiff[1] = grad_lambda[0][1] - grad_lambda[2][1];
				break;
			case 2:
				product = lambda[0] * lambda[1];
				diff = lambda[1] - lambda[0];
				dProd[0] = lambda[1] * grad_lambda[0][0] + lambda[0] * grad_lambda[1][0];
				dProd[1] = lambda[1] * grad_lambda[0][1] + lambda[0] * grad_lambda[1][1];
				dDiff[0] = grad_lambda[1][0] - grad_lambda[0][0];
				dDiff[1] = grad_lambda[1][1] - grad_lambda[0][1];
				break;
			default:
				throw std::runtime_error("Invalid edge index");
		}

		// grad phi^e_k = (lambda_1 * lambda_2)' * ker_{k-2} +
		// 				(lambda_1 * lambda_2) * ker_{k-2}' * (lambda_2 - lambda_1)'
		grad[0] = dProd[0] * kernel(k-2, diff) + product * kernel_deriv(k-2, diff) * dDiff[0];
		grad[1] = dProd[1] * kernel(k-2, diff) + product * kernel_deriv(k-2, diff) * dDiff[1];
		return;
	}

	// bubble functions
	int bubble_i = i - noVertices - noEdgeFuncs;
	
	// find n_1, n_2 corresponding to index
	int n1 = 1, n2 = 1, counter = 0;
	bool found = false;
	for (n1=1; n1<p-1; n1++)
	{
		for (n2=1; n2<p-n1; n2++)
		{
			if (counter == bubble_i)
			{
				found = true;
				break;
			}
			counter++;
		}
		if (found) break;
	}

	// find lambda_1 * lambda_2 * lambda_3 term and deriv
	// (lambda_1 * lambda_2 * lambda_3)' = lambda_1' * lambda_2 * lambda_3 +
	//									lambda_1 * lambda_2' * lambda_3 +
	// 									lambda_1 * lambda_2 * lambda_3'
	double prod = lambda[0] * lambda[1] * lambda[2];
	double dProd[2] = {
		grad_lambda[0][0] * lambda[1] * lambda[2] + lambda[0] * grad_lambda[1][0] * lambda[2] + lambda[0] * lambda[1] * grad_lambda[2][0],
		grad_lambda[0][1] * lambda[1] * lambda[2] + lambda[0] * grad_lambda[1][1] * lambda[2] + lambda[0] * lambda[1] * grad_lambda[2][1]
	};

	// find kernel terms and derivs
	double ker1 = kernel(n1-1, lambda[2] - lambda[1]);
	double dKer1 = kernel_deriv(n1-1, lambda[2] - lambda[1]);
	double ker2 = kernel(n2-1, lambda[1] - lambda[0]);
	double dKer2 = kernel_deriv(n2-1, lambda[1] - lambda[0]);

	// find deriv of difference terms in kernel
	double dDiff1[2] = {
		grad_lambda[2][0] - grad_lambda[1][0],
		grad_lambda[2][1] - grad_lambda[1][1]
	};
	double dDiff2[2] = {
		grad_lambda[1][0] - grad_lambda[0][0],
		grad_lambda[1][1] - grad_lambda[0][1]
	};

	// grad phi^b_{n1, n2} = (lambda_1 * lambda_2 * lambda_3)' * ker_{n1-1} * ker_{n2-1} +
	// 						(lambda_1 * lambda_2 * lambda_3) * ker_{n1-1}' * (lambda_3 - lambda_2)' * ker_{n2-1} +
	// 						(lambda_1 * lambda_2 * lambda_3) * ker_{n1-1} * ker_{n2-1}' * (lambda_2 - lambda_1)'
	grad[0] = dProd[0] * ker1 * ker2 +
				prod * dKer1 * ker2 * dDiff1[0] +
				prod * ker1 * dKer2 * dDiff2[0];
	grad[1] = dProd[1] * ker1 * ker2 +
				prod * dKer1 * ker2 * dDiff1[1] +
				prod * ker1 * dKer2 * dDiff2[1];
}