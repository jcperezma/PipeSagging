#include <vector>
#include "MathUtils.h"

using namespace std;


class StructuralFiniteElement {
public :
	virtual void  getElementArea(double &area, vector<Point2D> const & coords)=0;
	virtual void  computeLocalStifMatrix(double const &E, double const & nu,  vector<Point2D> const & coords)=0;
	virtual void  addLocalstifToGlobal(Matrix2D<double> & K )=0;
	void  setNodeIds(vector<int> & indices ) { 	nodes.assign(indices.begin(),indices.end());}
	void setThickness(double & thick) {this->thick = thick;}

protected:
	vector<int> nodes;		// Vector with node indices
	dMatrix2D<double> k;		// Element stifness matrix
	double thick;			// Element thickness
};



class linearTri : public StructuralFiniteElement {
public :
	//
	linearTri()  {}
	~linearTri() {}

	// overloaded functions
	void  getElementArea(double &area, vector<Point2D> const & coords);
	void  computeLocalStifMatrix(double const &E, double const & nu,  vector<Point2D> const & coords);
	void  addLocalstifToGlobal(Matrix2D<double> & K );

};


class linearQuad : public StructuralFiniteElement {
public :
	//
	linearQuad() {}
	~linearQuad() {}

	//overloaded functions
	void  getElementArea(double &area, vector<Point2D> const & coords);
	void  computeLocalStifMatrix(double const &E, double const & nu,  vector<Point2D> const & coords);
	void  addLocalstifToGlobal(Matrix2D<double> & K );

private :
	Matrix2D<double> EvaluateLocalK(const double &E, const double &nu, const double & Xi, const double &eta, const Point2D & p1, const Point2D & p2, const Point2D &  p3, const Point2D & p4  );
};