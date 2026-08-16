#include "GaussQuadrature.hpp"
#include "Polynomial.hpp"
#include <stdexcept>
#include <cmath>

// constructor
GaussQuadrature::GaussQuadrature()
 : points(), weights()
{
}

// destructor
GaussQuadrature::~GaussQuadrature()
{
}

void GaussQuadrature::assembleQuadrature(int n)
{
	// validate input
	if (n < 1)
	{
		throw std::invalid_argument("Number of quadrature points must be at least 1");
	}

	// precompute low order quadrature points and weights
	switch(n)
	{
		case 1:
			points = {0.0};
			weights = {2.0};
			return;
		case 2:
			points = {-1/sqrt(3), 1/sqrt(3)};
			weights = {1.0, 1.0};
			return;
		case 3:
			points = {-sqrt(3.0/5), 0.0, sqrt(3.0/5)};
			weights = {5.0/9, 8.0/9, 5.0/9};
			return;
		case 4:
			points = {-sqrt(3.0/7 + (2.0/7)*sqrt(6.0/5)), -sqrt(3.0/7 - (2.0/7)*sqrt(6.0/5)), sqrt(3.0/7 - (2.0/7)*sqrt(6.0/5)), sqrt(3.0/7 + (2.0/7)*sqrt(6.0/5))};
			weights = {(18.0 - sqrt(30))/36, (18.0 + sqrt(30))/36, (18.0 + sqrt(30))/36, (18.0 - sqrt(30))/36};
			return;
		case 5:
			// nodes are (1/3) sqrt(5 -+ 2 sqrt(10/7)), i.e. sqrt(5/9 -+ (2/9) sqrt(10/7))
			points = {-sqrt(5.0/9 + (2.0/9)*sqrt(10.0/7)), -sqrt(5.0/9 - (2.0/9)*sqrt(10.0/7)), 0.0, sqrt(5.0/9 - (2.0/9)*sqrt(10.0/7)), sqrt(5.0/9 + (2.0/9)*sqrt(10.0/7))};
			weights = {(322.0 - 13*sqrt(70))/900, (322.0 + 13*sqrt(70))/900, 128.0/225, (322.0 + 13*sqrt(70))/900, (322.0 - 13*sqrt(70))/900};
			return;
	}

	// for n > 5, use Newton's method to find roots of Legendre polynomial

	// resize vectors
	points.resize(n);
	weights.resize(n);

	// constants for Newton's method
	const double tol = 1e-14;
	const int max_iter = 100;

	// find points and weights
	for (int i=0; i<(n+1) / 2; i++)
	{
		// initial values for Newton's method
		double x = cos(M_PI * (i + 0.75) / (n + 0.5)); // initial guess using Chebyshev approximation
		double x_old;
		int iter = 0;

		// Newton's method main loop
		do
		{
			x_old = x;
			x = x_old - legendre(n, x_old) / legendre_deriv(n, x_old);
			iter++;
			if (iter >= max_iter)
			{
				// throw error if Newton's method did not converge
				throw std::runtime_error("Newton's method did not converge");
			}
		// fabs, not abs: the integer abs would truncate every step below 1.0 to zero
		// and stop the iteration after a single step
		} while (fabs(x - x_old) > tol);

		// store points and weights using symmetry of roots
		double w = 2.0 / ((1.0 - x * x) * pow(legendre_deriv(n, x), 2));
		points[i] = -x;
		weights[i] = w;
		if (i != n - i - 1)
		{
			points[n - i - 1] = x;
			weights[n - i - 1] = w;
		}
	}
}
