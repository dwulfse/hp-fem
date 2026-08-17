#include "FE_Solution.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// -------------------------------------------------------------------------
// benchmark problems
//
// each forcing function is paired with the analytic solution of the
// corresponding boundary value problem, so that the error can be measured
// -------------------------------------------------------------------------

// -u'' = pi^2 sin(pi x) on (0, 1), u(0) = u(1) = 0
double forcingTrig1D(const std::vector<double>& x)
{
	return M_PI * M_PI * sin(M_PI * x[0]);
}

double solutionTrig1D(const std::vector<double>& x)
{
	return sin(M_PI * x[0]);
}

// -u'' = 1 on (0, 1), u(0) = u(1) = 0
double forcingConst1D(const std::vector<double>&)
{
	return 1.0;
}

double solutionConst1D(const std::vector<double>& x)
{
	return x[0] * (1.0 - x[0]) / 2.0;
}

// -lap u = 2 pi^2 sin(pi x) sin(pi y), zero on the boundary
double forcingTrig2D(const std::vector<double>& x)
{
	return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

double solutionTrig2D(const std::vector<double>& x)
{
	return sin(M_PI * x[0]) * sin(M_PI * x[1]);
}

// -lap u = 1, zero on the boundary; no closed form on a general domain
double forcingConst2D(const std::vector<double>&)
{
	return 1.0;
}

// -------------------------------------------------------------------------
// configuration
// -------------------------------------------------------------------------

struct Config
{
	int dimension = 1;
	int degree = 1;
	int elements = 4;			// 1D only; the 2D element count comes from the mesh
	std::string mesh = "6";		// 2D only; loads domain.<mesh>.node and .ele
	std::string problem = "trig";
	int gridPoints = 64;
	bool field = false;
	int fieldLevels = 2;

	bool sweep = false;
	int sweepMaxDegree = 4;

	bool semilinear = false;
	int q = 1;
	double alpha = 0.7;
	double lambda = 1.0;
	int maxIter = 100;
	double tol = 1e-8;
};

void printUsage()
{
	std::cout <<
		"usage: FEM [options]\n"
		"\n"
		"solves -lap u = f, or -lap u + lambda u^(2q+1) = f, with zero dirichlet\n"
		"boundary conditions. must be run from main/, since mesh files are read\n"
		"and csv output written relative to the working directory.\n"
		"\n"
		"  -d, --dim <1|2>        spatial dimension (default 1)\n"
		"  -p, --degree <int>     polynomial degree (default 1)\n"
		"  -n, --elements <int>   number of elements, 1D only (default 4)\n"
		"  -m, --mesh <name>      2D mesh, loads domain.<name>.node (default 6)\n"
		"                         meshes 1-9 refine the unit square; L.1 is an\n"
		"                         L shaped domain\n"
		"      --problem <name>   trig or const (default trig)\n"
		"      --grid <int>       points per side written to solution.csv (default 64)\n"
		"      --field            also write field.csv, the triangulation and the\n"
		"                         solution on it, for plotting (2D only)\n"
		"      --field-levels <n> subdivisions per element in field.csv (default 2)\n"
		"\n"
		"      --sweep            also run an h-refinement study to hp_error.csv\n"
		"      --max-degree <int> highest degree in the sweep (default 4)\n"
		"\n"
		"      --semilinear       solve the semilinear problem instead (2D only)\n"
		"      --q <int>          exponent parameter, u^(2q+1) (default 1)\n"
		"      --alpha <float>    picard relaxation, 0 < alpha < 1 (default 0.7)\n"
		"      --lambda <float>   coefficient of the nonlinear term (default 1)\n"
		"      --max-iter <int>   picard iteration cap (default 100)\n"
		"      --tol <float>      picard convergence tolerance (default 1e-8)\n"
		"\n"
		"  -h, --help             show this message\n"
		"\n"
		"examples:\n"
		"  FEM                            1D, degree 1, 4 elements\n"
		"  FEM -d 2 -p 3 -m 6             2D on a 128 element mesh at degree 3\n"
		"  FEM -d 2 -p 2 -m L.1           2D on the L shaped domain\n"
		"  FEM -d 2 -p 3 --sweep          2D, plus the convergence study\n"
		"  FEM -d 2 --semilinear --q 2    semilinear with u^5\n";
}

namespace
{
	bool matches(const char* arg, const char* shortOpt, const char* longOpt)
	{
		return (shortOpt && strcmp(arg, shortOpt) == 0) || strcmp(arg, longOpt) == 0;
	}

	// reads the value following a flag, failing if it is missing
	const char* value(int argc, char** argv, int& i, const char* flag)
	{
		if (i + 1 >= argc)
		{
			std::cerr << "error: " << flag << " needs a value\n";
			exit(1);
		}
		return argv[++i];
	}
}

// returns false if the program should stop without solving, as for --help
bool parseArgs(int argc, char** argv, Config& cfg)
{
	for (int i = 1; i < argc; i++)
	{
		char* a = argv[i];

		if (matches(a, "-h", "--help"))
		{
			printUsage();
			return false;
		}
		else if (matches(a, "-d", "--dim"))
		{
			cfg.dimension = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, "-p", "--degree"))
		{
			cfg.degree = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, "-n", "--elements"))
		{
			cfg.elements = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, "-m", "--mesh"))
		{
			cfg.mesh = value(argc, argv, i, a);
		}
		else if (matches(a, nullptr, "--problem"))
		{
			cfg.problem = value(argc, argv, i, a);
		}
		else if (matches(a, nullptr, "--grid"))
		{
			cfg.gridPoints = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--field"))
		{
			cfg.field = true;
		}
		else if (matches(a, nullptr, "--field-levels"))
		{
			cfg.fieldLevels = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--sweep"))
		{
			cfg.sweep = true;
		}
		else if (matches(a, nullptr, "--max-degree"))
		{
			cfg.sweepMaxDegree = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--semilinear"))
		{
			cfg.semilinear = true;
		}
		else if (matches(a, nullptr, "--q"))
		{
			cfg.q = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--alpha"))
		{
			cfg.alpha = atof(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--lambda"))
		{
			cfg.lambda = atof(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--max-iter"))
		{
			cfg.maxIter = atoi(value(argc, argv, i, a));
		}
		else if (matches(a, nullptr, "--tol"))
		{
			cfg.tol = atof(value(argc, argv, i, a));
		}
		else
		{
			std::cerr << "error: unrecognised option " << a << "\n"
				<< "try FEM --help\n";
			exit(1);
		}
	}

	return true;
}

// rejects combinations the solver cannot handle, rather than failing obscurely later
bool validate(const Config& cfg)
{
	if (cfg.dimension != 1 && cfg.dimension != 2)
	{
		std::cerr << "error: dimension must be 1 or 2\n";
		return false;
	}
	if (cfg.degree < 1)
	{
		std::cerr << "error: polynomial degree must be at least 1\n";
		return false;
	}
	if (cfg.dimension == 1 && cfg.elements < 1)
	{
		std::cerr << "error: element count must be at least 1\n";
		return false;
	}
	if (cfg.problem != "trig" && cfg.problem != "const")
	{
		std::cerr << "error: problem must be trig or const\n";
		return false;
	}
	if (cfg.field && cfg.dimension != 2)
	{
		std::cerr << "error: field output is two dimensional only\n";
		return false;
	}
	if (cfg.semilinear && cfg.dimension != 2)
	{
		std::cerr << "error: the semilinear solver is two dimensional only\n";
		return false;
	}
	if (cfg.semilinear && (cfg.alpha <= 0.0 || cfg.alpha >= 1.0))
	{
		std::cerr << "error: the relaxation parameter must satisfy 0 < alpha < 1\n";
		return false;
	}

	return true;
}

// -------------------------------------------------------------------------
// drivers
// -------------------------------------------------------------------------

void runLinear(const Config& cfg)
{
	bool trig = (cfg.problem == "trig");

	double (*f)(const std::vector<double>&) =
		cfg.dimension == 1
			? (trig ? forcingTrig1D : forcingConst1D)
			: (trig ? forcingTrig2D : forcingConst2D);

	// constant forcing on a general 2D domain has no closed form, so there is
	// nothing to measure the error against in that one case
	double (*exact)(const std::vector<double>&) =
		cfg.dimension == 1
			? (trig ? solutionTrig1D : solutionConst1D)
			: (trig ? solutionTrig2D : nullptr);

	int n = (cfg.dimension == 1) ? cfg.elements : 1;
	FE_Solution FEM(n, cfg.degree, cfg.dimension);

	std::string mesh = (cfg.dimension == 2) ? cfg.mesh : "";
	std::vector<double> solution = FEM.solve(f, mesh);

	std::cout << "dimension  " << cfg.dimension << "\n"
		<< "degree     " << cfg.degree << "\n";
	if (cfg.dimension == 1)
	{
		std::cout << "elements   " << cfg.elements << "\n";
	}
	else
	{
		std::cout << "mesh       domain." << cfg.mesh << "\n";
	}
	std::cout << "mesh size  " << FEM.getMeshSize() << "\n"
		<< "dofs       " << solution.size() << "\n";

	if (exact)
	{
		std::cout << "L2 error   " << FEM.getL2Error(exact) << "\n";
		FEM.sendSolutionToFile(cfg.gridPoints, exact);
	}
	else
	{
		FEM.sendSolutionToFile(cfg.gridPoints);
	}

	if (cfg.field && cfg.dimension == 2)
	{
		FEM.sendFieldToFile(cfg.fieldLevels);
		std::cout << "wrote field.csv\n";
	}

	std::cout << "wrote solution.csv\n";
}

void runSemilinear(const Config& cfg)
{
	FE_Solution FEM(1, cfg.degree, 2);

	std::vector<double> solution = FEM.solveSemilinear(
		forcingTrig2D, cfg.q, cfg.alpha, cfg.lambda, cfg.maxIter, cfg.tol, cfg.mesh);

	std::cout << "degree     " << cfg.degree << "\n"
		<< "mesh       domain." << cfg.mesh << "\n"
		<< "exponent   u^" << (2 * cfg.q + 1) << "\n"
		<< "alpha      " << cfg.alpha << "\n"
		<< "lambda     " << cfg.lambda << "\n"
		<< "dofs       " << solution.size() << "\n";

	FEM.sendSolutionToFile(cfg.gridPoints);

	if (cfg.field)
	{
		FEM.sendFieldToFile(cfg.fieldLevels);
		std::cout << "wrote field.csv\n";
	}

	std::cout << "wrote solution.csv\n";
}

// h-refinement study: hold the degree fixed, refine the mesh, record the error
void runSweep(const Config& cfg)
{
	bool trig = (cfg.problem == "trig");

	std::ofstream file("hp_error.csv");
	file << "p,h,error\n";

	for (int p = 1; p <= cfg.sweepMaxDegree; p++)
	{
		if (cfg.dimension == 1)
		{
			double (*f)(const std::vector<double>&) = trig ? forcingTrig1D : forcingConst1D;
			double (*exact)(const std::vector<double>&) = trig ? solutionTrig1D : solutionConst1D;

			for (int n = 10; n <= 20; n++)
			{
				FE_Solution FEM(n, p, 1);
				FEM.solve(f);
				file << p << "," << FEM.getMeshSize() << "," << FEM.getL2Error(exact) << "\n";
			}
		}
		else
		{
			// domain.9 is excluded: it holds 15565 elements against domain.8's
			// 512, so it is a graded mesh rather than a uniform refinement of
			// the same family and does not belong in an h-refinement study
			const std::vector<std::string> meshes = {"4", "5", "6", "7", "8"};

			for (size_t i = 0; i < meshes.size(); i++)
			{
				FE_Solution FEM(1, p, 2);
				FEM.solve(forcingTrig2D, meshes[i]);
				double h = FEM.getMeshSize();
				double error = FEM.getL2Error(solutionTrig2D);
				file << p << "," << h << "," << error << "\n";
			}
		}

		std::cout << "swept degree " << p << "\n";
	}

	file.close();
	std::cout << "wrote hp_error.csv\n";
}

int main(int argc, char** argv)
{
	Config cfg;

	if (!parseArgs(argc, argv, cfg))
	{
		return 0;
	}
	if (!validate(cfg))
	{
		return 1;
	}

	if (cfg.semilinear)
	{
		runSemilinear(cfg);
	}
	else
	{
		runLinear(cfg);
	}

	if (cfg.sweep)
	{
		runSweep(cfg);
	}

	return 0;
}
