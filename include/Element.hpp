#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include <vector>
#include "PolynomialSpace.hpp"

class Element
{
	public:
		int id; // element id
		int p; // polynomial degree
		std::vector<int> local_DoF; // DoF
		PolynomialSpace poly;	// polynomial space

		// constructors
		Element() {};
		Element(int p, std::vector<int>& local_DoF) : p(p), local_DoF(local_DoF), poly(p) {};
		Element(int id, int p, std::vector<int>& local_DoF) : id(id), p(p), local_DoF(local_DoF), poly(p) {};

		// destructor
		virtual ~Element() {};

		// pure virtual methods
		virtual std::vector<std::vector<double>> getLocalStiffness() = 0;
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) = 0;
};

#endif