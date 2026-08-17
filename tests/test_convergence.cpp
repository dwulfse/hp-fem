// End to end tests: assemble and solve real problems, and check the answer.
//
// The tests above check individual pieces in isolation. These run the whole
// pipeline against problems with known analytic solutions and assert on the
// rate at which the error shrinks as the mesh is refined. Theory says the L2
// error behaves like O(h^(p+1)), so the slope of log(error) against log(h)
// should come out near p+1.
//
// A rate test is a strong end to end check precisely because it is a statement
// about a trend rather than a single number: a constant factor error in the
// assembly, the quadrature or the solve leaves the rate intact, but anything
// that breaks the approximation itself shows up immediately. Every one of the
// four defects this project carried moved one of these numbers.
//
// These tests read mesh files by relative path, so they run with main/ as the
// working directory, which is what the CMake setup arranges.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "FE_Solution.hpp"

#include <cmath>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
	// -u'' = pi^2 sin(pi x), with u(0) = u(1) = 0
	double forcingTrig1D(const std::vector<double>& x)
	{
		return M_PI * M_PI * sin(M_PI * x[0]);
	}

	double solutionTrig1D(const std::vector<double>& x)
	{
		return sin(M_PI * x[0]);
	}

	// -u'' = 1, with u(0) = u(1) = 0
	double forcingConst1D(const std::vector<double>&)
	{
		return 1.0;
	}

	double solutionConst1D(const std::vector<double>& x)
	{
		return x[0] * (1.0 - x[0]) / 2.0;
	}

	// -laplacian u = 2 pi^2 sin(pi x) sin(pi y), zero on the boundary
	double forcingTrig2D(const std::vector<double>& x)
	{
		return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]);
	}

	double solutionTrig2D(const std::vector<double>& x)
	{
		return sin(M_PI * x[0]) * sin(M_PI * x[1]);
	}

	// least squares slope of log(error) against log(h)
	double fittedRate(const std::vector<double>& h, const std::vector<double>& error)
	{
		size_t n = h.size();
		double mx = 0.0, my = 0.0;
		for (size_t i = 0; i < n; i++)
		{
			mx += log(h[i]) / n;
			my += log(error[i]) / n;
		}

		double num = 0.0, den = 0.0;
		for (size_t i = 0; i < n; i++)
		{
			num += (log(h[i]) - mx) * (log(error[i]) - my);
			den += (log(h[i]) - mx) * (log(h[i]) - mx);
		}
		return num / den;
	}
}

TEST_CASE("1D solution matches the analytic solution pointwise", "[solver][1d]")
{
	// four elements at degree three is enough to be visually exact
	FE_Solution FEM(4, 3, 1);
	FEM.solve(forcingTrig1D);

	for (int i = 0; i <= 4; i++)
	{
		double x = i / 4.0;
		INFO("x = " << x);
		REQUIRE_THAT(FEM.evaluateSolution({x}), WithinAbs(solutionTrig1D({x}), 5e-4));
	}
}

TEST_CASE("1D reproduces a quadratic solution to machine precision", "[solver][1d]")
{
	// u = x(1-x)/2 is a quadratic, so it lies in the p >= 2 space exactly and
	// the finite element solution should be the analytic one up to round off
	FE_Solution FEM(4, 2, 1);
	FEM.solve(forcingConst1D);

	for (int i = 0; i <= 8; i++)
	{
		double x = i / 8.0;
		INFO("x = " << x);
		REQUIRE_THAT(FEM.evaluateSolution({x}), WithinAbs(solutionConst1D({x}), 1e-12));
	}
}

TEST_CASE("1D attains the theoretical convergence rate", "[solver][1d][convergence]")
{
	for (int p = 1; p <= 4; p++)
	{
		std::vector<double> h, error;

		for (int n = 10; n <= 20; n++)
		{
			FE_Solution FEM(n, p, 1);
			FEM.solve(forcingTrig1D);
			h.push_back(FEM.getMeshSize());
			error.push_back(FEM.getL2Error(solutionTrig1D));
		}

		double rate = fittedRate(h, error);
		INFO("p = " << p << ", fitted rate " << rate << ", expected " << p + 1);
		REQUIRE_THAT(rate, WithinAbs(p + 1.0, 0.15));
	}
}

TEST_CASE("2D attains the theoretical convergence rate", "[solver][2d][convergence]")
{
	// domain.9 is deliberately absent from this list: it is a graded mesh
	// rather than a uniform refinement of the others, so it does not belong in
	// an h-refinement study
	const std::vector<std::string> meshes = {"4", "5", "6", "7", "8"};

	for (int p = 1; p <= 4; p++)
	{
		std::vector<double> h, error;

		for (const std::string& mesh : meshes)
		{
			FE_Solution FEM(1, p, 2);
			FEM.solve(forcingTrig2D, mesh);
			h.push_back(FEM.getMeshSize());
			error.push_back(FEM.getL2Error(solutionTrig2D));
		}

		double rate = fittedRate(h, error);
		INFO("p = " << p << ", fitted rate " << rate << ", expected " << p + 1);
		REQUIRE_THAT(rate, WithinAbs(p + 1.0, 0.2));
	}
}

TEST_CASE("refining the mesh reduces the error monotonically", "[solver][2d][convergence]")
{
	// weaker than a rate test but it localises a regression: if this passes and
	// the rate test fails, the approximation is still improving, just too slowly
	const std::vector<std::string> meshes = {"4", "5", "6", "7", "8"};

	for (int p = 1; p <= 3; p++)
	{
		double previous = -1.0;

		for (const std::string& mesh : meshes)
		{
			FE_Solution FEM(1, p, 2);
			FEM.solve(forcingTrig2D, mesh);
			double error = FEM.getL2Error(solutionTrig2D);

			if (previous > 0.0)
			{
				INFO("p = " << p << ", mesh " << mesh);
				REQUIRE(error < previous);
			}
			previous = error;
		}
	}
}

TEST_CASE("the mesh size of a uniform 1D mesh is 1/n", "[mesh][1d]")
{
	for (int n = 4; n <= 16; n *= 2)
	{
		FE_Solution FEM(n, 1, 1);
		FEM.solve(forcingConst1D);
		INFO("n = " << n);
		REQUIRE_THAT(FEM.getMeshSize(), WithinRel(1.0 / n, 1e-12));
	}
}

TEST_CASE("2D mesh sizes shrink by root two between successive refinements", "[mesh][2d]")
{
	// each mesh in the family doubles the element count, so h falls by sqrt(2)
	// rather than by 2; assuming otherwise silently halves every measured rate
	const std::vector<std::string> meshes = {"4", "5", "6", "7", "8"};
	std::vector<double> h;

	for (const std::string& mesh : meshes)
	{
		FE_Solution FEM(1, 1, 2);
		FEM.solve(forcingTrig2D, mesh);
		h.push_back(FEM.getMeshSize());
	}

	for (size_t i = 0; i + 1 < h.size(); i++)
	{
		INFO("meshes " << meshes[i] << " -> " << meshes[i + 1]);
		REQUIRE_THAT(h[i] / h[i + 1], WithinRel(sqrt(2.0), 1e-6));
	}
}

TEST_CASE("the semilinear solver converges", "[solver][semilinear]")
{
	// no analytic solution is available here, so this asserts only that the
	// Picard iteration terminates and returns something finite and bounded
	FE_Solution FEM(1, 2, 2);
	std::vector<double> u = FEM.solveSemilinear(forcingTrig2D, 1, 0.7, 1.0, 100, 1e-8, "6");

	REQUIRE(u.size() > 0);

	double largest = 0.0;
	for (size_t i = 0; i < u.size(); i++)
	{
		INFO("dof " << i);
		REQUIRE(std::isfinite(u[i]));
		largest = std::max(largest, fabs(u[i]));
	}

	// the forcing is the same as the linear problem whose solution peaks at 1,
	// and adding a positive reaction term can only damp it
	REQUIRE(largest > 0.1);
	REQUIRE(largest < 2.0);
}
