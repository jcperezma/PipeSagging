#include "FiniteElement.h"

double computeVisc(double const & T, double const & zero_shear_mu, double const & alpha){
	double visc =0.0;
 		 
	if(T >120 ){
		visc= zero_shear_mu * exp(alpha/T  );
		return visc;
	}

	if(T<=120 && T >100){
		double x2 = 120;
		double y2 = zero_shear_mu * exp(alpha/120  );	// I could pre compute this
		double x1 = 100;
		double y1 = y2*100;
		double m = (y2-y1)/(x2-x1);
		double b = y1 - m * x1;
		visc = m * T +b;
		return visc;
	}
	if(T<= 100){
		visc = zero_shear_mu * exp(alpha/120  )*100;
	}

	return visc;
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

void computeNiandDel(double const & xi, double const & eta, dMatrix2D<double> & Del, vector<double>& Ni){
	// computes shape functions and derivatives of shpae functions evaluated and xi and eta. 
	// there are n*n number of Ni and Del arrays, can be computed apriori. I would be saving ~50 operationn per integration
	// point, which means 50xn^2 operations fewer per element, 450 for 9 noded elements, 
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
	Ni.resize(9);
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

	Del.resize(2,9);
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

}

vector<vector<integrationPointVals>> initializeIntegrationVector(){
	
	vector<vector<integrationPointVals>> gaussPointVals(3);
	// for the cases where you use 1, 2 or 3 gauss points
	for (int i = 1; i < 4; i++) 
	{
		vector <double>   wGauss;
		vector <double>   gaussPoint;
		initGaussVars(   i,  wGauss, gaussPoint);

		for (int m = 0; m < i; m++)
		{
			for (int n = 0; n < i; n++)
			{
				// compute Ni, Del and w1w2 for this combination of gauss points
				integrationPointVals val;
				val.w1w2 = wGauss[m] * wGauss[n]; 
				computeNiandDel(gaussPoint[m],gaussPoint[n],val.del,val.Ni);
				gaussPointVals[i-1].push_back(val);					
			}
		}

	}
	return gaussPointVals;
}


vector<vector<integrationPointVals>> BiLinearQuadThermalFluid::gaussPointVals = initializeIntegrationVector();


void BiLinearQuadThermalFluid::computeDetandDerivs(vector<Point2D> const & coords, dMatrix2D<double>   Del, vector<double>const  & Ni, double &detJ, dMatrix2D<double>&B){
	
	// Compute Jacobian matrix
	// J = |  dx/dxi  dy/dxi  |
	//	   |  dx/deta dy/deta |	

	dMatrix2D<double>  J(2,2);     // Jacobian Matrix
	dMatrix2D<double>  J_inv(2,2); // Inverse of the Jacobian
	J = 0;

	for (int i = 0; i < 9; i++)
	{
		int NodeID = nodes[i]; // global node ID
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
	 B =  J_inv.matMul(Del);
	
}


void BiLinearQuadThermalFluid::getIntegrationData(int ngauss, int gaussPointID,double & w1w2, dMatrix2D<double> & Del, vector<double> & Ni){ 
		w1w2 = gaussPointVals[ngauss-1][gaussPointID].w1w2;  
		Del = gaussPointVals[ngauss-1][gaussPointID].del; 
		Ni =gaussPointVals[ngauss-1][gaussPointID].Ni; };

void BiLinearQuadThermalFluid::computeFluidK(vector<Point2D> const & coords,vector<double> & U, vector<double>& T, dMatrix2D<double> & k, vector<double> & f, SimParameters & params){
	// For fluids there are two degrees of freedom per node, 
	// the solution u is organized as u = { u1, u2,..un, v1, v2,..vn }, n is the number of nodes i.e. coords.size()
	// the load vector is organized as f = {fx1,fx2,..fxn, fy1, fy2,..fyn}
	// the local stiffnes matrix k is 18x18
	// the local load vector f is 18x1



	// the indexing here might be confusing, 
	// Pittman breaks k into for smaller matrices 9x9
	// k = | A11 A12 |
	//	   | A21 A12 |
	// and the force vector f into
	// f_x f_y
	// Ill first do it filling the four submatrices and subvectors following my notation
	//

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
	dMatrix2D<double>  del;
	double detJ;

	// external parameters, I duplicate the values here to not make the code below more verbose.
	//double mu = params.mu;
	double lambda = params.lambda*params.mu;
	double rho =params.rho;
	double g_y =params.g_y;
	//double mu=43;

	// find viscosity at center point

	double T_i = T[nodes[8]];
    	double mu = computeVisc(T_i,params.mu,params.alpha);

	
	for (int i = 0; i < nnodes; i++) 
	{
		f[i*2] = 0; 
		f[i*2+1] = 0; 
	}

	// integrate with respect of xi and eta
	int pointID=-1;
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			pointID++;
			// find shape functions and derivatives at xi and eta
			
			//this->computeShapeFunctions(gaussPoint[xi],gaussPoint[eta], coords, B, Ni,  detJ);
			
			double w1w2;
			//vector<double> Ni2;
			getIntegrationData(nGauss,pointID,w1w2,del,Ni);
			computeDetandDerivs(coords,del,Ni,detJ,B);
			double detWeight = detJ * w1w2; //Combined weights and detJ

			// find viscosity at the integration point, i.e interpolate temperature 
			/*
			double T_i = 0;
				for (int i=0;i<nnodes;i++) 
				{
					double T_par = T[nodes[i]];
					T_i+= T[nodes[i]]*Ni[i];
				}
				

				double mu = params.mu;*/
				// mu = computeVisc(T_i,params.mu,params.alpha);
			//double detWeight = detJ * wGauss[xi] * wGauss[eta]; //Combined weights and detJ
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
	nGauss--;
	gaussPoint.clear();wGauss.clear();
	initGaussVars( nGauss, wGauss, gaussPoint);
	dMatrix2D<double> k_lambda(nnodes*2,nnodes*2);
	k_lambda =0;

	pointID=-1;
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			pointID++;
			// find shape functions and derivatives at xi and eta
			double w1w2;
			//vector<double> Ni2;
			getIntegrationData(nGauss,pointID,w1w2,del,Ni);
			computeDetandDerivs(coords,del,Ni,detJ,B);
			double detWeight = detJ * w1w2; //Combined weights and detJ
			// in my notation the i and j are switched, in pittman's paper it would be swapped
			double T_i = 0;
				for (int i=1;i<nnodes;i++) T_i+= T[nodes[i]]*Ni[i];

				//double mu = params.mu;
				//lambda = params.lambda* computeVisc(T_i,params.mu,params.alpha);

			for (int i = 0; i < nnodes; i++)
			{
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

void BiLinearQuadThermalFluid::computeThermalK(vector<Point2D> const & coords,vector<double> & U, vector<double>& T, dMatrix2D<double> & k, vector<double> & f, SimParameters & params){
	// the local stiffnes matrix k is 18x18
	// the local load vector f is 18x18

	dMatrix2D<double> k_visc(nnodes,nnodes);
	vector<double>f_x(nnodes); // RHS
	k=0;

	// We have to integrate the values over the element, we will use gauss quadrature with three points in each direction
	int nGauss = params.nGauss;
	// initialize weight and variable vectors

	vector<double> wGauss, gaussPoint;
	initGaussVars( nGauss, wGauss, gaussPoint);

	dMatrix2D<double>  B;
	dMatrix2D<double>  del;
	vector<double>  Ni(nnodes);
	double detJ;

	// external parameters, I duplicate the values here to not make the code below more verbose.
	double kappa = params.kappa;

	
	for (int i = 0; i < nnodes; i++) 
	{
		f_x[i] = 0; 
	}
	int pointID =-1;
	// integrate with respect of xi and eta
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			pointID++;
			// find shape functions and derivatives at xi and eta
			double w1w2;
			//vector<double> Ni2;
			getIntegrationData(nGauss,pointID,w1w2,del,Ni);
			computeDetandDerivs(coords,del,Ni,detJ,B);
			double detWeight = detJ * w1w2; //Combined weights and detJ
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
				}
			}
		}
	}


}


void BiLinearQuadThermalFluid::computeThermalKandM(vector<Point2D> const & coords,vector<double> & U, vector<double>& T, dMatrix2D<double> & k,dMatrix2D<double> & m, vector<double> & f, SimParameters & params){
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
	dMatrix2D<double>  del;
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

	k=0;
	m=0;

	int pointID = -1;
	// integrate with respect of xi and eta
	for (int xi= 0; xi < nGauss; xi++)
	{
		for (int eta = 0; eta < nGauss; eta++)
		{
			pointID++;
			// find shape functions and derivatives at xi and eta
			double w1w2;
			//vector<double> Ni2;
			getIntegrationData(nGauss,pointID,w1w2,del,Ni);
			computeDetandDerivs(coords,del,Ni,detJ,B);
			double detWeight = detJ * w1w2; //Combined weights and detJ
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
void BiLinearQuadThermalFluid::computeLocalStifMatrix(vector<Point2D> const & coords, vector<double> & U, vector<double>& T, dMatrix2D<double> & k, vector<double> & f, physType const & type, SimParameters & params){

	switch (type)
	{
	case FLUID :
		computeFluidK(coords,U,T, k, f, params);
			  break;
	case THERMAL :
		computeThermalK(coords,U,T, k, f, params);
			  break;
	default:
		break;
	}


}

void BiLinearQuadThermalFluid::computeLocalStifandMassMatrix(vector<Point2D> const & coords,vector<double> & U, vector<double>& T, dMatrix2D<double> & k , dMatrix2D<double> & m, vector<double> & f, physType const & type, SimParameters & params){

	switch (type)
	{
	case FLUID :
		//computeFluidKandM(coords, k, m, f, params);
			  break;
	case THERMAL :
		computeThermalKandM(coords,U,T, k,m, f, params);
			  break;
	default:
		break;
	}


}


void  BiLinearQuadThermalFluid::addLocalstifToGlobal( dMatrix2D<double>  & k, vector<double> const & f, FE_data & data, physType const & type){
	// For FLUID system k is 18x18
	// For THERMAL system k is 9x9

	for (int i = 0; i < nnodes; i++)
	{
		for (int m = 0; m < data.dim; m++)
		{
			data.F[this->nodes[i]*data.dim+m] += f[i*data.dim+m];

			for (int n = 0; n < data.dim; n++)
			{
				for (int j = 0; j < nnodes; j++)
				{
					int a = data.dim * this->nodes[i]+m;  int b =data.dim * this->nodes[j]+n;

					data.K[a][b] += k[data.dim *i+m][data.dim*j+n];
				}
			}

		}
	}

}


void  BiLinearQuadThermalFluid::addLocalStifandMassToGlobal( dMatrix2D<double>  & k,dMatrix2D<double> & mass, vector<double> const & f, FE_data & data, physType const & type){
	// For FLUID system k is 18x18
	// For THERMAL system k is 9x9

	for (int i = 0; i < nnodes; i++)
	{
		for (int m = 0; m < data.dim; m++)
		{
			data.F[this->nodes[i]*data.dim+m] += f[i*data.dim+m];

			for (int n = 0; n < data.dim; n++)
			{
				for (int j = 0; j < nnodes; j++)
				{
					int a = data.dim * this->nodes[i]+m;  int b =data.dim * this->nodes[j]+n;

					data.K[a][b] += k[data.dim *i+m][data.dim*j+n];
					data.M[a][b] += mass[data.dim*i+m][data.dim*j+n];
				}
			}

		}
	}

}