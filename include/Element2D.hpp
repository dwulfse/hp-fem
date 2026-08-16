#ifndef ELEMENT2DHEADERDEF
#define ELEMENT2DHEADERDEF

#include <vector>
#include "Element.hpp"
#include "Point2D.hpp"

class Element2D : public Element
{
	public:
		std::vector<Point2D> nodes;

		// constructors
		Element2D();
		Element2D(int id, int p, std::vector<int>& local_DoF, std::vector<Point2D>& nodes);

		// destructor
		virtual ~Element2D();

		// methods
		double basisSign(int i);
		virtual std::vector<std::vector<double>> getLocalStiffness() override;
		std::vector<double> getLocalStiffnessProduct(const std::vector<double>& U);
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) override;
		std::vector<double> getLocalNonlinearLoad(const std::vector<double>& U, int q);
};

#endif