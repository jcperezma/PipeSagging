#include "FiniteElement.h"



void BiLinearQuadThermalFluid::computeShapeFunctions(double const & xi, double const & eta, vector<Point2D> const & coords, dMatrix2D<double> & B, vector<double>& Ni, double & detJ){
	// computes shape functions and derivatives of shpae functions evaluated and xi and eta. 

	// computes:
	// matrix B
	// B =   | dNi/dx | 0
	//	     | dNi/dy | 1
	// Shape functions Ni(xi,eta) 
	// Determinant of the Jacobian
	
	double minuseta  = 0.5*eta*(1-eta);
	double pluseta   = 0.5*eta*(1+eta);
	double minusxi   = 0.5*xi*(1-xi);
	double plusxi    = 0.5*xi*(1+xi);
	double minuseta2 = 1-eta*eta;
	double minusxi2  = 1-xi*xi;
	
	// compute shape functions Ni(xi,eta)

	Ni[0]=  minusxi  *  minuseta;
	Ni[1]= -plusxi   *  minuseta;
	Ni[2]=  plusxi   *  pluseta;
	Ni[3]= -minusxi  *  pluseta;
	Ni[4]= -minusxi2 *  minuseta;
	Ni[5]=  plusxi   *  minuseta2;
	Ni[6]=  minusxi2 *  pluseta;
	Ni[7]= -minusxi  *  minuseta2;
	Ni[8]=  minusxi2 *  minuseta2;

	// compute derivatives of shape functions
	// Del = | dNi/dxi  | 0
	//	     | dNi/deta | 1
	
	double minus2xi  = 0.5*(1-2*xi);
	double plus2xi   = 0.5*(1+2*xi);
	double minus2eta = 0.5*(1-2*eta);
	double plus2eta  = 0.5*(1+2*eta);

	dMatrix2D<double> Del(2,9);
	// dNi/dxi
	Del[0][0] =  minus2xi * minuseta ;
	Del[0][1] = -plus2xi  * minuseta;
	Del[0][2] =  plus2xi  * pluseta;
	Del[0][3] = -minus2xi * pluseta;
	Del[0][4] =  2* xi    * minuseta;
	Del[0][5] =  plus2xi  * minuseta2;
	Del[0][6] = -2*xi     * pluseta;
	Del[0][7] = -minus2xi * minuseta2;
	Del[0][8] =  -2*xi     * minuseta2;

	// dNi/deta
	Del[1][0] =  minusxi   * minus2eta ;
	Del[1][1] = -plusxi    * minus2eta;
	Del[1][2] =  plusxi    * plus2eta;
	Del[1][3] = -minusxi   * plus2eta;
	Del[1][4] = -minusxi2  * minus2eta;
	Del[1][5] = -2*plusxi  * eta; 
	Del[1][6] =  minusxi2  * plus2eta;
	Del[1][7] =  2*minusxi * eta;
	Del[1][8] = - 2*minusxi2* eta;

	// Compute Jacobian matrix
	// J = |  dx/dxi  dy/dxi  |
	//	   |  dx/deta dy/deta |	

	dMatrix2D<double>  J(2,2);     // Jacobian Matrix
	dMatrix2D<double>  J_inv(2,2); // Inverse of the Jacobian
	J = 0;

	for (int i = 0; i < 9; i++)
	{
		int NodeID = this->nodes[i]; // global node ID
		Point2D curPosition = coords[NodeID];
		double partial = Del[0][i]*coords[NodeID].x;
		J[0][0] += Del[0][i]*coords[NodeID].x ; // dx/dxi  = sum (dNi/dxi Xi)
		J[0][1] += Del[0][i]*coords[NodeID].y ; // dy/dxi  = sum (dNi/dxi Yi)
		J[1][0] += Del[1][i]*coords[NodeID].x ; // dx/deta = sum (dNi/deta Xi)
		J[1][1] += Del[1][i]*coords[NodeID].y ; // dy/deta = sum (dNi/deta Yi)
	}

	detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0] ;

	// Compute inverse of the Jacobian
	double coeff = 1 / detJ;
	J_inv[0][0] = coeff *  J[1][1];
	J_inv[0][1] = coeff * -J[0][1];
	J_inv[1][0] = coeff * -J[1][0];
	J_inv[1][1] = coeff *  J[0][0];

	// compute matrix B
	// B =   | dNi/dx | 0
	//	     | dNi/dy | 1

	B = J_inv.matMul(Del);

}

void initGaussVars(int const &   nGauss, vector<double> & wGauss,vector <double> &  gaussPoint){
	// initialize the gauss variables used for the integration depending on how many gauss points will be used
	switch (nGauss)
	{
	case 1:
		wGauss.push_back(2); 
		gaussPoint.push_back(0); 
		break;
	case 2:
		wGauss.push_back(1); wGauss.push_back(1); 
		gaussPoint.push_back(0.57735); gaussPoint.push_back(-0.57735); 
		break;
	case 3:
		wGauss.push_back(0.5555); wGauss.push_back(0.5555); wGauss.push_back(0.8888);
		gaussPoint.push_back(-0.774596); gaussPoint.push_back(0.774596); gaussPoint.push_back(0);
		break;
	default:
		break;
	}
}

void BiLinearQuadThermalFluid::computeFluidKandM(vector<Point2D> const & coords, vector<double>const & T, vector<double>const & U,  dMatrix2D<double> & k, dMatrix2D<double> & m,vector<double> & f, SimParameters & params, bool computeM){
	// For fluids there are two degrees of freedom per node, 
	// the solution u is organized as u = { u1, u2,..un, v1, v2,..vn }, n is the number of nodes i.e. coords.size()
	// the load vector is organized as f = {fx1,fx2,..fxn, fy1, fy2,..fyn}
	// the local stiffnes matrix k is 18x18
	// the local load vector f is 18x1

		
	dMatrix2D<double> k_visc(nnodes*2,nnodes*2);
	vector<double>f_x(nnodes);
	vector<double>f_y(nnodes);


	// We have to integrate the values over the element, we will use gauss quadrature with three points in each direction
	int nGauss = params.nGauss;
	// initialize weight and variable vectors

	vector<double> wGauss, gaussPoint;
	initGaussVars( nGauss, wGauss, gaussPoint);


	
	dMatrix2D<double>  B;
	vector<double>  Ni(nnodes);
	double detJ;

	// external parameters, I duplicate the values here to not make the code below more verbose.
	double mu = params.mu;
	double lambda = params.lambda;
	double rho =params.rho;
	double g_y =params.g_y;

	
	for (int i = 0; i < nnodes; i++) 
	{
		f[i*2] = 0; 
		f[i*2+1] = 0; 
	}

	// integrate with respect of xi and eta
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			// find shape functions and derivatives at xi and eta
			this->computeShapeFunctions(gaussPoint[xi],gaussPoint[eta], coords, B, Ni,  detJ);
			double detWeight = detJ * wGauss[xi] * wGauss[eta]; //Combined weights and detJ
			// in my notation the i and j are switched, in pittman's paper it would be swapped
			for (int j = 0; j < nnodes; j++)
			{
				for (int i = 0; i < nnodes; i++)
				{
					// remember
					// B =   | dNi/dx | 0
					//	     | dNi/dy | 1

					// with stokes 
					enum {dx, dy};
					k_visc[2*i][2*j] += mu*(2 * B[dx][i] * B[dx][j] + B[dy][i] * B[dy][j])* detWeight;
					k_visc[2*i+1][2*j+1] += mu*(B[dx][i] * B[dx][j] + 2 * B[dy][i] * B[dy][j])* detWeight;
					k_visc[2*i][2*j+1] += mu*( B[dy][i] * B[dx][j])* detWeight;
					k_visc[2*i+1][2*j] += mu*( B[dx][i] * B[dy][j])* detWeight;

				}
			}
		}
	}

	// find penalty matrix, with a reduced integration
	nGauss = 2;
	gaussPoint.clear();wGauss.clear();
	initGaussVars( nGauss, wGauss, gaussPoint);
	dMatrix2D<double> k_lambda(nnodes*2,nnodes*2);
	k_lambda =0;

	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			// find shape functions and derivatives at xi and eta
			this->computeShapeFunctions(gaussPoint[xi],gaussPoint[eta], coords, B, Ni,  detJ);
			double detWeight = detJ * wGauss[xi] * wGauss[eta]; //Combined weights and detJ
			// in my notation the i and j are switched, in pittman's paper it would be swapped
			for (int i = 0; i < nnodes; i++)
			{
				f[i*2] = 0;							  // there are no source terms in X
				f[i*2+1] += -rho * g_y *Ni[i]*detWeight; // source term in Y is gravity
				for (int j = 0; j < nnodes; j++)
				{
					// remember
					// B =   | dNi/dx | 0
					//	     | dNi/dy | 1
					enum {dx, dy};

					k_lambda[2*i][2*j] += lambda*(B[dx][i] * B[dx][j] )* detWeight;
					k_lambda[2*i+1][2*j+1] += lambda*(B[dy][i] * B[dy][j])* detWeight;
					
					k_lambda[2*i][2*j+1] += lambda*( B[dx][i] * B[dy][j])* detWeight;
					k_lambda[2*i+1][2*j] += lambda*( B[dy][i] * B[dx][j])* detWeight;

				}
			}
		}
	}

	k = k_visc + k_lambda;
			
}


void BiLinearQuadThermalFluid::computeThermalKandM(vector<Point2D> const & coords, vector<double>const & T, vector<double>const & U, dMatrix2D<double> & k,dMatrix2D<double> & m, vector<double> & f, SimParameters & params, bool computeM){
	// the local stiffnes matrix k is 18x18
	// the local load vector f is 18x18

	dMatrix2D<double> k_visc(nnodes,nnodes);
	vector<double>f_x(nnodes); // RHS


	// We have to integrate the values over the element, we will use gauss quadrature with three points in each direction
	int nGauss = params.nGauss;
	// initialize weight and variable vectors

	vector<double> wGauss, gaussPoint;
	initGaussVars( nGauss, wGauss, gaussPoint);

	dMatrix2D<double>  B;
	vector<double>  Ni(nnodes);
	double detJ;

	// external parameters, I duplicate the values here to not make the code below more verbose.
	double kappa = params.kappa;
	double Cp =params.Cp;
	double rho = params.rho;
	
	for (int i = 0; i < nnodes; i++) 
	{
		f_x[i] = 0; 
	}

	// integrate with respect of xi and eta
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			// find shape functions and derivatives at xi and eta
			this->computeShapeFunctions(gaussPoint[xi],gaussPoint[eta], coords, B, Ni,  detJ);
			double detWeight = detJ * wGauss[xi] * wGauss[eta]; //Combined weights and detJ
			// in my notation the i and j are switched, in pittman's paper it would be swapped
			for (int i = 0; i < nnodes; i++)
			{
				f_x[i] = 0;							  // there are no source terms 

				for (int j = 0; j < nnodes; j++)
				{
					// remember
					// B =   | dNi/dx | 0
					//	     | dNi/dy | 1
					// without convection there is only difussion, and it is a linear constitutive relationship
					enum {dx, dy};
					k[i][j] += kappa * (B[dx][i]*B[dx][j]+B[dy][i]*B[dy][j]) *detWeight; 

					//mass matrix
					m[i][j] += rho * Cp * ( Ni[i]*Ni[j]) * detWeight;
				}
			}
		}
	}


}


void BiLinearQuadThermalFluid::computeLocalStifandMassMatrix(vector<Point2D> const & coords,vector<double> const & U,vector<double> const & T, dMatrix2D<double> & k, dMatrix2D<double> & m,  vector<double> & f, physType const & type, SimParameters & params, bool computeM){

	switch (type)
	{
	case FLUID :
		computeFluidKandM(coords,U,T, k,m, f, params, computeM);
			  break;
	case THERMAL :
		computeThermalKandM(coords,U,T, k,m, f, params, computeM);
			  break;
	default:
		break;
	}


}


void  BiLinearQuadThermalFluid::addLocalStifandMassToGlobal( dMatrix2D<double>  & k,dMatrix2D<double> & m, vector<double> const & f, FE_data & data, physType const & type, bool computeM){
	// For FLUID system k is 18x18
	// For THERMAL system k is 9x9

	int totalnodes; // total nodes in the simulation
	
	totalnodes = data.F.size()/data.dim;


	for (int i = 0; i < nnodes; i++)
	{
		data.F[this->nodes[i]] += f[i*data.dim];

		if (type == FLUID) data.F[this->nodes[i]+totalnodes] += f[i*data.dim+1];

		for (int j = 0; j < nnodes; j++)
		{
			// If the system is fluid need to add more elements to K
			if (type == FLUID)
			{
				int rx = 2 * this->nodes[i];  int sx =2 * this->nodes[j];
				int ry = 2 * this->nodes[i]+1; int sy =2 * this->nodes[j]+1;
				data.K[rx][sx] += k[2*i][2*j];
				data.K[rx][sy] += k[2*i][2*j+1];
				data.K[ry][sx] += k[2*i+1][2*j];
				data.K[ry][sy] += k[2*i+1][2*j+1];
				// TODO addd mass to fluid, dont need it for pipe sagging
			}
			if (type == THERMAL)
			{
				data.K[this->nodes[i]][this->nodes[j]] += k[i][j];
				if(computeM)
					data.M[this->nodes[i]][this->nodes[j]] += m[i][j];
			}

		}
	}

}