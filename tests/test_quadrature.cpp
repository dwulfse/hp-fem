// Gaussian quadrature.
//
// The defining property of an n-point Gauss-Legendre rule is that it integrates
// every polynomial up to degree 2n-1 exactly. That single statement is enough to
// pin down the nodes and weights completely, which makes it an unusually good
// thing to test: any error in either produces a failure here.
//
// Both quadrature bugs this project carried would have been caught by the first
// test below. The 5-point node table had one wrong digit, and the Newton
// iteration used the integer abs() so it stopped after a single step.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "GaussQuadrature.hpp"
#include "GaussQuadrature2D.hpp"

#include <cmath>

using Catch::Matchers::WithinAbs;

namespace
{
	// exact value of the integral of x^k over [-1, 1]
	double exactInterval(int k)
	{
		return (k % 2 == 1) ? 0.0 : 2.0 / (k + 1);
	}

	// exact value of the integral of xi^a eta^b over the reference triangle
	// {-1 < xi, eta ; xi + eta < 0}
	double exactTriangle(int a, int b)
	{
		double s = ((b + 1) % 2 == 0 ? 1.0 : -1.0) / (b + 1);
		return s * (exactInterval(a + b + 1) - exactInterval(a));
	}
}

TEST_CASE("1D Gauss rules are exact to degree 2n-1", "[quadrature]")
{
	for (int n = 1; n <= 10; n++)
	{
		GaussQuadrature quad;
		quad.assembleQuadrature(n);

		REQUIRE(static_cast<int>(quad.points.size()) == n);
		REQUIRE(static_cast<int>(quad.weights.size()) == n);

		for (int k = 0; k <= 2 * n - 1; k++)
		{
			double approx = 0.0;
			for (int i = 0; i < n; i++)
			{
				approx += quad.weights[i] * pow(quad.points[i], k);
			}

			INFO("n = " << n << ", integrating x^" << k);
			REQUIRE_THAT(approx, WithinAbs(exactInterval(k), 1e-12));
		}
	}
}

TEST_CASE("1D Gauss nodes lie inside the interval and weights are positive", "[quadrature]")
{
	// a rule with a node outside [-1, 1] or a negative weight is not a Gauss
	// rule, and would evaluate basis functions outside the reference element
	for (int n = 1; n <= 10; n++)
	{
		GaussQuadrature quad;
		quad.assembleQuadrature(n);

		double total = 0.0;
		for (int i = 0; i < n; i++)
		{
			INFO("n = " << n << ", node " << i);
			REQUIRE(quad.points[i] > -1.0);
			REQUIRE(quad.points[i] < 1.0);
			REQUIRE(quad.weights[i] > 0.0);
			total += quad.weights[i];
		}

		// the weights must sum to the length of the interval
		REQUIRE_THAT(total, WithinAbs(2.0, 1e-12));
	}
}

TEST_CASE("triangle rules integrate low order monomials exactly", "[quadrature]")
{
	// the triangle rule is a collapsed tensor product of the 1D rule, so it
	// inherits the 1D accuracy but the Duffy factor costs one degree
	for (int nq = 2; nq <= 8; nq++)
	{
		GaussQuadrature2D quad;
		quad.assembleQuadrature(nq);

		REQUIRE(quad.points.size() == quad.weights.size());
		REQUIRE(quad.points.size() > 0);

		for (int deg = 0; deg <= nq - 1; deg++)
		{
			for (int a = 0; a <= deg; a++)
			{
				int b = deg - a;
				double approx = 0.0;
				for (size_t i = 0; i < quad.points.size(); i++)
				{
					approx += quad.weights[i]
						* pow(quad.points[i].x, a)
						* pow(quad.points[i].y, b);
				}

				INFO("nq = " << nq << ", integrating xi^" << a << " eta^" << b);
				REQUIRE_THAT(approx, WithinAbs(exactTriangle(a, b), 1e-11));
			}
		}
	}
}

TEST_CASE("triangle rule weights sum to the reference triangle area", "[quadrature]")
{
	// the reference triangle {-1 < xi, eta ; xi + eta < 0} has area 2
	for (int nq = 2; nq <= 8; nq++)
	{
		GaussQuadrature2D quad;
		quad.assembleQuadrature(nq);

		double total = 0.0;
		for (size_t i = 0; i < quad.weights.size(); i++)
		{
			total += quad.weights[i];
		}

		INFO("nq = " << nq);
		REQUIRE_THAT(total, WithinAbs(2.0, 1e-12));
	}
}
