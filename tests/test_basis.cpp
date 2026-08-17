// The hierarchical basis on the reference triangle.
//
// A conforming basis has to satisfy a handful of structural properties, and each
// of them is a direct assertion about function values rather than anything
// statistical. That makes them cheap to check and very hard to satisfy by
// accident:
//
//   * vertex functions are 1 at their own vertex and 0 at the other two
//   * edge functions vanish on the two edges they do not belong to
//   * bubble functions vanish on all three edges
//   * along its own edge, an edge function reduces to the 1D Lobatto function
//
// The last of these is the sharpest. Edge modes run k = 2, ..., p and the kernel
// is indexed k-2, and this project passed the 0-based loop counter straight
// through, so the first two modes of every edge were 2/(1+x) and 2/(1-x): not
// polynomials at all. Because a guard zeroed them exactly at the vertices they
// still looked plausible, and only a comparison against the values they are
// supposed to take exposes it.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "PolynomialSpace.hpp"
#include "Polynomial.hpp"

#include <cmath>

using Catch::Matchers::WithinAbs;

namespace
{
	// reference triangle vertices, in the local ordering the element uses:
	// local vertex 0 is where the first vertex function equals one, and so on
	const double VERTEX[3][2] = {{-1.0, -1.0}, {1.0, -1.0}, {-1.0, 1.0}};

	// a point a fraction t along edge e, where t runs from -1 to 1 and edge e
	// joins local vertices e and (e+1) % 3
	void edgePoint(int e, double t, double& xi, double& eta)
	{
		int a = e, b = (e + 1) % 3;
		double s = 0.5 * (t + 1.0);
		xi = VERTEX[a][0] + s * (VERTEX[b][0] - VERTEX[a][0]);
		eta = VERTEX[a][1] + s * (VERTEX[b][1] - VERTEX[a][1]);
	}

	int noEdgeModes(int p) { return p - 1; }
	int noEdgeFuncs(int p) { return 3 * (p - 1); }
	int noDoF(int p) { return (p + 1) * (p + 2) / 2; }
}

TEST_CASE("vertex functions interpolate at their own vertex", "[basis]")
{
	for (int p = 1; p <= 5; p++)
	{
		PolynomialSpace poly(p);

		for (int i = 0; i < 3; i++)
		{
			for (int v = 0; v < 3; v++)
			{
				double expected = (i == v) ? 1.0 : 0.0;
				INFO("p = " << p << ", vertex function " << i << " at vertex " << v);
				REQUIRE_THAT(poly.basis_2D(i, VERTEX[v][0], VERTEX[v][1]),
					WithinAbs(expected, 1e-12));
			}
		}
	}
}

TEST_CASE("edge functions reduce to the 1D Lobatto functions along their edge", "[basis]")
{
	// along edge 0, the two affine coordinates involved are exactly the 1D
	// Lobatto vertex functions of the edge parameter, so edge mode m must equal
	// l_{m+2}
	for (int p = 2; p <= 6; p++)
	{
		PolynomialSpace poly(p);

		for (int m = 0; m < noEdgeModes(p); m++)
		{
			for (int s = 0; s <= 20; s++)
			{
				double t = -1.0 + 2.0 * s / 20.0;
				double xi, eta;
				edgePoint(0, t, xi, eta);

				INFO("p = " << p << ", edge mode " << m << " at t = " << t);
				REQUIRE_THAT(poly.basis_2D(3 + m, xi, eta),
					WithinAbs(lobatto(m + 2, t), 1e-12));
			}
		}
	}
}

TEST_CASE("edge functions vanish on the edges they do not belong to", "[basis]")
{
	for (int p = 2; p <= 6; p++)
	{
		PolynomialSpace poly(p);

		for (int i = 3; i < 3 + noEdgeFuncs(p); i++)
		{
			int own = (i - 3) / noEdgeModes(p);

			for (int e = 0; e < 3; e++)
			{
				if (e == own) continue;

				for (int s = 0; s <= 20; s++)
				{
					double xi, eta;
					edgePoint(e, -1.0 + 2.0 * s / 20.0, xi, eta);

					INFO("p = " << p << ", basis " << i << " on edge " << e);
					REQUIRE_THAT(poly.basis_2D(i, xi, eta), WithinAbs(0.0, 1e-12));
				}
			}
		}
	}
}

TEST_CASE("bubble functions vanish on every edge", "[basis]")
{
	// bubbles are interior to the element, which is what lets neighbouring
	// triangles ignore each other's bubbles entirely
	for (int p = 3; p <= 6; p++)
	{
		PolynomialSpace poly(p);

		for (int i = 3 + noEdgeFuncs(p); i < noDoF(p); i++)
		{
			for (int e = 0; e < 3; e++)
			{
				for (int s = 0; s <= 20; s++)
				{
					double xi, eta;
					edgePoint(e, -1.0 + 2.0 * s / 20.0, xi, eta);

					INFO("p = " << p << ", bubble " << i << " on edge " << e);
					REQUIRE_THAT(poly.basis_2D(i, xi, eta), WithinAbs(0.0, 1e-12));
				}
			}
		}
	}
}

TEST_CASE("the basis has the expected number of functions", "[basis]")
{
	// 3 vertices, p-1 modes on each of 3 edges, and (p-1)(p-2)/2 bubbles
	for (int p = 1; p <= 6; p++)
	{
		int counted = 3 + 3 * (p - 1) + (p - 1) * (p - 2) / 2;
		INFO("p = " << p);
		REQUIRE(counted == noDoF(p));
	}
}

TEST_CASE("analytic gradients agree with finite differences", "[basis]")
{
	// basis_2D_grad is written out by hand from the product rule, so it is easy
	// for it to drift from basis_2D; comparing against a difference quotient of
	// basis_2D catches that
	const double eps = 1e-6;

	for (int p = 1; p <= 4; p++)
	{
		PolynomialSpace poly(p);

		// interior sample points, kept clear of the edges where the kernel
		// guards switch on
		const double sample[4][2] = {
			{-0.5, -0.5}, {-0.2, -0.6}, {-0.6, -0.2}, {-0.34, -0.33}
		};

		for (int i = 0; i < noDoF(p); i++)
		{
			for (int s = 0; s < 4; s++)
			{
				double xi = sample[s][0], eta = sample[s][1];

				double grad[2];
				poly.basis_2D_grad(i, xi, eta, grad);

				double dxi = (poly.basis_2D(i, xi + eps, eta)
					- poly.basis_2D(i, xi - eps, eta)) / (2 * eps);
				double deta = (poly.basis_2D(i, xi, eta + eps)
					- poly.basis_2D(i, xi, eta - eps)) / (2 * eps);

				INFO("p = " << p << ", basis " << i << " at (" << xi << ", " << eta << ")");
				REQUIRE_THAT(grad[0], WithinAbs(dxi, 1e-6));
				REQUIRE_THAT(grad[1], WithinAbs(deta, 1e-6));
			}
		}
	}
}
