#include <vector>
#include "MathUtils.h"
using namespace std;

enum physType{FLUID, THERMAL};

// Finite Element Data
struct FE_data{
	vector<double> U;		// Displacements array, if 2D array (u1, u2, ..., un, v1, v2, ...vn)
	vector<double> U_old;	// Old displacements
	vector<int> DOFid;		// holds the indices for all the velocities degrees of freedom
	vector<int> BC_F_id;	// index of force applied to mesh at node, the index is given in the global DOF vector notation, i.e. F = {Fx1, Fy1, Fx2, Fy2, ...,FxN, FyN}, same for all the other arrays
	vector<double> BC_F;	// Values of force applied at nodes
	vector<int> BC_U_id;	// Physical constrain indices
	vector<double> BC_U;	// Values of physical constrains
	spMatrix2D<double> K;	
	spMatrix2D<double> K_reduced; // must be same type as K
	spMatrix2D<double> M;	// mass Matrix
	vector<double>F;
	double U_o;				// initial displacement field	
	int dim;				// num of Dimensions of variable U
};




// simulation parameters needed by the thermal fluid element
struct SimParameters{
	double rho;     // fluid density
	double mu;	    // fluid viscosity
	double lambda;  // penalty constant
	double g_y;     // gravity
	int nGauss;     // number of Gauss points used in the integration
	double kappa;	// thermal conductivity
	double Cp;		// heat capacity
	double dt;		//time step
	double theta;	// theta parameter

};


struct integrationPointVals{
	dMatrix2D<double> del; // derivatives of the shape functions in local coords
	vector<double> Ni;		// shape functions
	double w1w2;			// combined gauss point weight
};

class FiniteElement {
public :
	virtual void  getElementArea(double &area, vector<Point2D> const & coords)=0;
	virtual void  computeLocalStifMatrix(vector<Point2D> const & coords, dMatrix2D<double> & k, vector<double> & f, physType const & type, SimParameters & params)=0;
	virtual void  computeLocalStifandMassMatrix(vector<Point2D> const & coords, dMatrix2D<double> & k, dMatrix2D<double> & m,  vector<double> & f, physType const & type, SimParameters & params)=0;
	virtual void  addLocalstifToGlobal( dMatrix2D<double> & k, vector<double> const & f, FE_data & data, physType const & type)=0;
	virtual void  addLocalStifandMassToGlobal( dMatrix2D<double> & k,dMatrix2D<double> & m, vector<double> const & f, FE_data & data, physType const & type)=0;

	void  setNodeIds(vector<int> & indices ) { 	nodes.assign(indices.begin(),indices.end());}
	vector<int> getNodes(){return nodes;};

protected:
	vector<int> nodes;		// Vector with node indices
	// no need to store k or f
	//Matrix2D<double> k;		// Element stifness matrix
	//vector<double> f;		// Vector with loads
};



class BiLinearQuadThermalFluid : public FiniteElement {
	/* Node winding
4-------7-------3
|		         |
|				 |
8       9 		 6
|				 |
|				 |
1-------5-------2
 */
public:
	 BiLinearQuadThermalFluid() {nnodes =9;};
	~BiLinearQuadThermalFluid() {}
	void  getElementArea(double &area, vector<Point2D> const & coords){};
	virtual void  computeLocalStifMatrix(vector<Point2D> const & coords, dMatrix2D<double> & k, vector<double> & f, physType const & type, SimParameters & params);
	virtual void  computeLocalStifandMassMatrix(vector<Point2D> const & coords, dMatrix2D<double> & k, dMatrix2D<double> & m,  vector<double> & f, physType const & type, SimParameters & params);
	virtual void  addLocalstifToGlobal( dMatrix2D<double> & k, vector<double> const & f, FE_data & data, physType const & type) override;
	virtual void  addLocalStifandMassToGlobal( dMatrix2D<double> & k,dMatrix2D<double> & m, vector<double> const & f, FE_data & data, physType const & type) override;

	//void computeShapeFunctions(double const & Xi, double const & eta, vector<Point2D> const & coords, dMatrix2D<double> & B, vector<double>& Ni,double & detJ);
	
	//static void intitializeIntegrationData();

protected:
	void computeFluidK(vector<Point2D> const & coords, dMatrix2D<double> & k, vector<double> & f, SimParameters & params);
	void computeThermalK(vector<Point2D> const & coords, dMatrix2D<double> & k, vector<double> & f, SimParameters & params);
	//void computeFluidKandM(vector<Point2D> const & coords, dMatrix2D<double> & k, dMatrix2D<double> & m,vector<double> & f, SimParameters & params);
	void computeThermalKandM(vector<Point2D> const & coords, dMatrix2D<double> & k, dMatrix2D<double> & m,vector<double> & f, SimParameters & params);
	int nnodes ; //number of nodes per element
	void getIntegrationData(int ngauss, int gaussPointID,double & w1w2, dMatrix2D<double> & Del, vector<double> & Ni);
	void computeDetandDerivs(vector<Point2D> const & coords, dMatrix2D<double>Del, vector<double>const & Ni, double &detJ,dMatrix2D<double>&B );
	static vector<vector<integrationPointVals>> gaussPointVals;
};

