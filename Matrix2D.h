#pragma once
#include <list>
#include <vector>
#include "Row.h"
#include <iostream>
//#include "MathUtils.h"
using namespace std;


template <typename T>
class Matrix2D{
	// Very simple Matrix class, allows to do simple operations such as matrix multiplication, addition, transpose
	// and solve using gauss seidel or SOR
public:
	// 
	//Matrix2D(){};
	virtual row<T>& operator[](int i)=0;
	virtual const row<T>& operator[](int i)const=0;
	virtual void resize(int nRows, int nColumns)=0;
	virtual void removeZeros() =0;
	virtual Matrix2D<T>&  transpose() =0 ;
	virtual void addRow(row<T> & row)= 0;
	~Matrix2D() {}
	
	virtual row<T> & getRow(int i) =0;
	vector<T>const operator*(vector<T>const & b){
		vector<T> c;
		for (int i = 0; i < getnRows(); i++)
		{
			T part= getRow(i)*b;
			c.push_back(part);
		}
		return c;

	};

	// Math operations
	/*
	virtual Matrix2D<T>& operator= ( const Matrix2D<T> &A)=0; // has to be implemented 
	virtual Matrix2D<T>& operator= (const T & value)=0; // set all elements to a certain value
	virtual	vector<T> dot(vector<T> b) =0; // pure virtual
	virtual Matrix2D<T> transpose()=0; //pure virtual
	virtual Matrix2D<T>& operator*= (const T & value)=0;
	virtual Matrix2D<T> operator* (const T & value)=0;
	virtual Matrix2D<T> operator * (const Matrix2D<T> & B)=0;
	virtual Matrix2D<T> operator + (const Matrix2D<T> & B)=0;
	virtual	Matrix2D<T> & operator += (const Matrix2D<T> & B)=0;
	*/
	/*
	virtual int getNRows()=0;
	virtual int getNColumns() const {return nColumns;}
	*(
	/*
	void resize(int newRow, int newColumn){nRows = newRow; nColumns= newColumn ;
		data.resize(nRows);
		for (int i = 0; i < nRows; i++)
		{
			data[i].resize(nColumns);
		}
	};*/



	
	virtual int const getnRows()const=0;
protected:
	Matrix2D(int n):nRows(n){};
	void setnRows(int newRows){nRows=newRows;};	
	int nRows;
};


// dense Matrix
template<typename T>
class dMatrix2D : public Matrix2D<T>{

public:
	dMatrix2D():Matrix2D(0){};

	dMatrix2D(int nRows, int nColumns ) : Matrix2D(nRows) , nColumns(nColumns) {
		matrixData.resize(nRows);
		for (int i = 0; i < nRows; i++)
		{
			matrixData[i].resize(nColumns);
		}
		
	};
	
	dMatrix2D<T>( dMatrix2D<T>const &A):Matrix2D(0){
	
		nRows = A.matrixData.size();
		nColumns = A.getNColumns();
		matrixData.resize(nRows);
		//resize matrix
		for (int i = 0; i < nRows; i++)
		{
			matrixData[i].resize(nColumns);
		}
		//Fill Matrix
		for (int i = 0; i < nRows; i++)
			{
				for (int j = 0; j < nColumns; j++)
				{
					matrixData[i][j] = A[i][j];
				}
			}
	}

	virtual void addRow(row<T> & row)override{
		
	};
	row<T>& operator[](int i){
		return matrixData[i];
	};

	const row<T>& operator[](int i)const override{
		return matrixData[i];
	};

	virtual void resize(int newRows, int newCols){
		matrixData.resize(newRows);
		for (int i = 0; i < newRows; i++)
		{
			matrixData[i].resize(newCols);
		}

		setnRows(newRows);
		nColumns = newCols;
	};

	virtual row<T> & getRow(int i){return matrixData[i]; };

	Matrix2D<T> &  transpose() {
		static dMatrix2D<T> trans(nColumns,nRows);
		for (int i = 0; i < nRows; i++)
			{
				for (int j = 0; j < nColumns; j++)
				{
					trans[j][i] = matrixData[i][j];
				}
			}
	return trans;
	}



	virtual int getNColumns() const {return nColumns;}


	dMatrix2D<T>& operator= (const T & value){ 
		for (int i =0; i< getnRows(); i++)
		{
			for (int j =0; j< getNColumns(); j++)
				{
					matrixData[i][j]=value;
				}
		}
	return *this;
	}

	/*
	virtual Matrix2D<T>& operator= (const T & value){};
	virtual Matrix2D<T>& operator= ( const Matrix2D<T> &A){} ;
	*/

	void removeZeros(){

	};
	dMatrix2D<T> matMul(dMatrix2D<T>const & B) {

	// Matrix multiplication
	dMatrix2D<T> C(nRows,B.getNColumns());
		C =0.0;
		if (nColumns == B.getnRows())
		{
			for (int i = 0; i < nRows; i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{

					for (int k = 0; k < nColumns; k++)
					{
						C[i][j] += matrixData[i][k]*B[k][j];
					}
				}
			}

		}
		return C;
	};

	/*
	virtual Matrix2D<T> transpose();
	virtual vector<T> dot(vector<T> b);

	virtual Matrix2D<T>& operator*= (const T & value);
	virtual	Matrix2D<T> operator* (const T & value);
	virtual Matrix2D<T> operator* (const Matrix2D<T> & B);
	virtual Matrix2D<T> operator + (const Matrix2D<T> & B);
	virtual Matrix2D<T> & operator += (const Matrix2D<T> & B);
	*/
	dMatrix2D<T> operator* (const T & value){ 
		// Scalar multiplication
		dMatrix2D<T> C(nRows,nColumns);
		C =0.0;
		
			for (int i = 0; i < nRows; i++)
			{
				for (int j = 0; j < nColumns; j++)
				{
					C[i][j] = operator[](i)[j]*value;
				}
			}		
		return C;
	}

	const dMatrix2D<T> & transpose_() {
	
		dMatrix2D<T> trans(getNColumns(),getnRows());
		for (int i = 0; i < getnRows(); i++)
			{
				for (int j = 0; j < getNColumns(); j++)
				{
					trans[j][i] = matrixData[i][j];
				}
			}
	return trans;
}
	void reduce(dMatrix2D<T> & K_reduced , vector<int> const & DOF_ID, vector<int> const & BC_ID){
		int num_displacement_DOF = DOF_ID.size();
		K_reduced.resize(num_displacement_DOF,num_displacement_DOF);
		K_reduced=0;
	for (int i = 0; i < num_displacement_DOF; i++)
		{
			for (int j = 0; j < num_displacement_DOF; j++)
			{
				K_reduced[i][j] = matrixData[DOF_ID[i]][DOF_ID[j]];
			}
		}
	}
	virtual int const getnRows()const override  { return matrixData.size();};

	protected:
		int nColumns;
		vector<Denserow<T>> matrixData;
};



// Sparse Matrix
template<typename T>

class spMatrix2D : public Matrix2D<T>{
		
public:
	spMatrix2D(int n=0): Matrix2D(n){
		resize(n,0); 
	}; // constructor

	
	row<T>& operator[](int i){
		return matrixData[i];
	};
	const row<T>& operator[](int i)const override{
		return matrixData[i];
	};

	void addItem(int const i, int const j, T const & entry){
		matrixData[i].insertNextItem(entry,j);
		
	};

	void addItem(int const i, int const j, T const & entry, bool isBC){
		matrixData[i].insertNextItem(entry,j, isBC);
		
	};

	void setToZero(){
		for (int i = 0; i < matrixData.size(); i++)
		{
			matrixData[i].setToZero();
		}
	};

	virtual void addRow(row<T> & row)override{
	};
	void addRow(Sparserow<T> & row){
		matrixData.push_back(row);
	};
	void removeZeros(){
		for (int i = 0; i < nRows; i++)
		{
			matrixData[i].removeZeros();
		}
	};

	void removeBC(vector<int> & idMap){
		for (int i = 0; i < matrixData.size(); i++)
		{
			matrixData[i].removeBCcolumns(idMap);
		}
	};
	/*
	spMatrix2D<T>& spMatrix2D<T>::operator= ( const spMatrix2D<T> &A){
	
		nRows = A.getnRows();
		matrixData.resize(nRows);
		//resize matrix
		//Fill Matrix
		for (int i = 0; i < nRows; i++)
			{
					matrixData[i] = A[i];
			}
		return *this;
	
	}*/

	void resize(int n, int m){
		
		matrixData.resize(n); 
	};

	virtual row<T> & getRow(int i){ return matrixData[i];  };
	~spMatrix2D(){};

	 Matrix2D<T> & transpose(){
		static spMatrix2D<T> tranposedMat(getnRows());
		for (int i = 0; i < getnRows(); i++)
		{
			for (int j = 0; j < matrixData[i].getSize(); j++)
			{
				rowElement<T> re = matrixData[i].getElementAt(j);
				if(re.getColumn()!=-1)
					tranposedMat.addItem(re.getColumn(),i,re.getValue());
			}
			
		}

	return tranposedMat;
}

 friend void reduce(spMatrix2D<T> & K, spMatrix2D<T> & K_reduced , vector<int> & DOF_ID, vector<int> & BC_ID)
{
		int num_BC = BC_ID.size();
		

		// indices to map the new column Ids
		vector<int> idMap(K.matrixData.size());

		int count =0;
		int indexCount=1;
		for (int i = 0; i < K.matrixData.size(); i++)
		{
			if (count<BC_ID.size()&&i==BC_ID[count])
			{
				
				count++; // get next BC_ID
			}
			idMap[i] = count;
		}

		K_reduced.matrixData.clear();
		K_reduced.resize(DOF_ID.size(),DOF_ID.size());

		for (int i =0; i <DOF_ID.size();  i++)
		{
			K_reduced.matrixData[i] = K.matrixData[DOF_ID[i]];
			//K_reduced.addRow(matrixData[DOF_ID[i]]);
		}

		K_reduced.removeBC(idMap);

		

	}
	 
	 void eraseRow(int rowID){
		 matrixData.erase(matrixData.begin()+rowID);
	 }

	

	 virtual int const getnRows()const override { return matrixData.size();};
protected :
	 vector<Sparserow<T>> matrixData;
};


template <typename T>
T  dotVectors(vector<T> const & a, vector<T> const & b  ){
	// assumes both vectors have same length and type
	T result = 0;

	for (int i = 0; i < a.size(); i++)
	{
		result += a[i] *b[i];
	}

	return result;
}


template <typename T>
vector<T> backSubs(Matrix2D<T> & U, vector<T> const & b  ){
	// ForwardSubstitution
	
	int n = b.size();
	vector<T> result(n);
	for (int i = n-1; i >=0 ; i--)
	{
		result[i] =(b[i]); //y[i] = b[i];
		
		result[i] -= (U[i].partialMultFrom(result, i+1));

			/*
		for (int j = i+1; j < n; j++)
		{
			result[i] -= U.data[i][j] * result[j]; //y[i] -= L[i, j] * y[j];
		}
		*/
		result[i] /= (U[i][i]); //y[i] /= L[i, i];
		
	}

	return result;
}

template <typename T>
void computeIcompleteChol( Matrix2D<T> &M){
	
	int n = M.getnRows();

	// for all rows
	for (int k = 0; k < n; k++)
	{
		// diagonal element 
		M[k][k]= sqrt(M[k][k]);

		// elements to the right of the diagonal
		for (int i = k+1; i < n; i++) // identical
		{
			if (M[i][k]!=0)  // if this exists  in the sparse matrix
				M[i][k]=M[i][k] /M[k][k];  // row operation          
		}

		// elements below the diagonal
		for (int j = k+1; j < n; j++) // terribly inneficient for the sparse
		{
			for (int i = j; i < n; i++)
			{
				if (M[i][j]!=0) // identical
					M[i][j] -= M[i][k]*M[j][k];
			}
		}


	}


	for (int i = 0; i < n; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			if (M[i][j]!=0)
				M[i][j] =0;
		}
	}


}

template <typename T>
vector<T> forwardSubs(Matrix2D<T>  & L, vector<T> const & b  ){
	// ForwardSubstitution
	vector<T> result;
	int n = b.size();
	for (int i = 0; i < n; i++)
	{
		result.push_back(b[i]); //y[i] = b[i];

		result[i] -= L[i].partialMultTo(result, i);
		/*
		for (int j = 0; j < i; j++)
		{
			result[i] -= L.data[i][j] * result[j]; //y[i] -= L[i, j] * y[j];
		}*/
		result[i] /= L[i][i]; //y[i] /= L[i, i];
	}

	return result;
}

// some useful Vector operations
template<typename T>
T operator * (vector<T> const & a, vector<T> const & b ){
	T result = 0;
	for (int i = 0; i < a.size(); i++)
		result+= a[i]*b[i];	 
	return result;
}

template<typename T>
vector<T> operator * (vector<T> const & a, T const & b ){
	vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i]= a[i]*b;	 
	return result;
}

template<typename T>
vector<T> operator * ( T const & b, vector<T> const & a ){
	vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i]*b;	 
	return result;
}

template <typename T >
vector<T> operator - (vector<T> const & a, vector<T> const & b ){
	vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i]-b[i];	 

	return result;
}

template <typename T >
vector<T> operator + (vector<T> const & a, vector<T> const & b ){
	vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i]+b[i];	 
	return result;
}




template<typename T>
vector<T> matMul (Matrix2D<T> & a, vector<T> const & b ){
	int n = a.getnRows();
	vector<T> result(n);
	for (int i = 0; i < n; i++)
	{
		result[i] = a[i] * b; // row times vector, each is implemented
	}
return result;			
}

template<typename T>
dMatrix2D<T> matMul (dMatrix2D<T> & A, dMatrix2D<T>  & B ){
	// Matrix multiplication
	Matrix2D<T> C(A.getnRows(),B.getNColumns());
		C =0.0;
		if (A.getNColumns() == B.getNRows())
		{
			for (int i = 0; i < A.getnRows(); i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{
					for (int k = 0; k < nColumns; k++)
					{
						C.data[i][j] += A[i][k]*B[k][j];
					}
				}
			}

		}
		return C;
	
}
template <typename T>
dMatrix2D<T> operator + (dMatrix2D<T> & A, dMatrix2D<T>  & B ){
	// Matrix multiplication
	dMatrix2D<T> C(A.getnRows(),A.getNColumns());
		C =0.0;
			for (int i = 0; i < A.getnRows(); i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{
						C[i][j] = A[i][j]+B[i][j];
				}
			}
		return C;
}
		
template <typename T>
dMatrix2D<T> operator - (dMatrix2D<T> & A, dMatrix2D<T>  & B ){
	// Matrix multiplication
	Matrix2D<T> C(A.getnRows(),A.getNColumns());
		C =0.0;
			for (int i = 0; i < A.getnRows(); i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{
						C.data[i][j] = A[i][j]-B[i][j];
				}
			}

		return C;
}

template <typename T>
dMatrix2D<T> operator * (dMatrix2D<T> & A, T B ){
	// Matrix multiplication
	Matrix2D<T> C(A.getnRows(),A.getNColumns());
		for (int i = 0; i < A.getnRows(); i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{
						C.data[i][j] = A[i][j]*B;
				}
			}

		return C;
}

template <typename T>
dMatrix2D<T> operator * ( T B,dMatrix2D<T> & A ){
	// Matrix multiplication
	Matrix2D<T> C(A.getnRows(),A.getNColumns());
		for (int i = 0; i < A.getnRows(); i++)
			{
				for (int j = 0; j < B.getNColumns(); j++)
				{
						C.data[i][j] = A[i][j]*B;
				}
			}

		return C;
}

template <typename T>
spMatrix2D<T> operator + (spMatrix2D<T> & A, spMatrix2D<T>  & B ){
	// 
	spMatrix2D<T> C=A;

			for (int i = 0; i < A.getnRows(); i++)
			{
				A[i].rowPlus(B[i],C[i]);
			}
		return C;
}
		
template <typename T>
spMatrix2D<T> operator - (spMatrix2D<T> & A, spMatrix2D<T>  & B ){
	// substraction
	// Assumes A and B have the same non zero elements 
	spMatrix2D<T> C=A;

			for (int i = 0; i < A.getnRows(); i++)
			{
				A[i].rowMinus(B[i],C[i]);
			}
		return C;
}

template <typename T>
spMatrix2D<T> operator * (spMatrix2D<T> & A, T B ){
	// substraction
	// Assumes A and B have the same non zero elements 
	spMatrix2D<T> C=A;

			for (int i = 0; i < A.getnRows(); i++)
			{
				A[i].rowTimes(B,C[i]);
			}

		return C;
}

template <typename T>
spMatrix2D<T> operator * ( T B,spMatrix2D<T> & A ){
	// substraction
	// Assumes A and B have the same non zero elements 
	spMatrix2D<T> C=A;

			for (int i = 0; i < A.getnRows(); i++)
			{
				A[i].rowTimes(B,C[i]);
			}

		return C;
}

// Solvers-------------------------------//


// Gauss Seidel and SOR
template <typename T>
	vector<T> solveSystemIter(Matrix2D<T>&A,vector<T> &b){
		// Uses Gauss Seidel to solve the linear system of equations AX=b, 
		// where A is the current matrix
		// assumes that b is the correct size

		// initialize guess
		vector<T> x;
		for (size_t  i = 0; i < b.size(); i++) 	x.push_back(0);

		int maxIter =100;
		double tol = 0.001;
		double error = 1000000000;
		int iter =1;
		int nRows = A.getnRows();
		while (iter<maxIter && error>tol)
		{
			// Update guesses
			for (int i = 0; i < nRows; i++)
			{
				double lambda = 1.7;  //if lambda > 1 over relaxation, lambda < 1 under relaxation 1 Gauss-Seidel
				//if (iter>25) {lambda =1;}
				x[i] = (1-lambda)*x[i] + (lambda/A[i][i]) * ( b[i] - (A[i]*x) + A[i][i]*x[i]  ); // SOR 
			}

			// compute result
			vector<T> result= A * x ;
		
			// compute error
			error =0;
			for (size_t  i = 0; i < result.size(); i++) error += abs(result[i]-b[i]);

			std::cout<<iter << " "<<error<<std::endl;
			iter++;
		}

		return x;
	}

	template <typename T>
	vector<T> solveSystemIterOld(Matrix2D<T> & A, vector<T> &b){
		// Uses Gauss Seidel to solve the linear system of equations AX=b, 
		// where A is the current matrix
		// assumes that b is the correct size

		// initialize guess
		vector<T> x;
		for (size_t  i = 0; i < b.size(); i++) 	x.push_back(0);

		int maxIter =A.getnRows();
		double tol = 0.001;
		double error = 10000000;
		int iter =1;
		while (iter<maxIter && error>tol)
		{
			// Update guesses
			for (int i = 0; i <A.getnRows(); i++)
			{
				T sum1= 0;
				for (int j = 0; j < A.getnRows(); j++) sum1+= A[i][j]*x[j];
			
				//x[i] = (1/data[i][i]) * ( b[i] - sum1 + data[i][i]*x[i]  ); // gauss-seidel
				
				double lambda = 1.3;
				if (iter>25) {lambda =1;}
				x[i] = (1-lambda)*x[i] + (lambda/A[i][i]) * ( b[i] - sum1 + A[i][i]*x[i]  ); // SOR
			}

			// compute result
			vector<T> result= A * x ;
		
			// compute error
			error =0;
			for (size_t  i = 0; i < result.size(); i++) error += abs(result[i]-b[i]);
			std::cout<<iter << " "<<error<<std::endl;
			iter++;
		}

		return x;
	}

	// Conjugate Gradient
	template <typename T >
vector<T> solveSystemIterCG(Matrix2D<T> & A, vector<T> &b){
	// uses conjugate gradient method 
	T tol = 0.1;
	// initialize guess
	vector<T> x;
	for (size_t  i = 0; i < b.size(); i++) 	x.push_back(0);

	// vector of residuals
	vector<T> r =b - (A*x); //r = b - A * x;
	
	// rsold = r' * r;
	T rsold=r*r;

	//p = r;
	vector<T> p = r;
	
	// iterate
	int iter =0;
	int maxIter = b.size();

	while (iter<maxIter && sqrt(rsold)>tol)
	{
		//Ap = A * p;
		vector<T> Ap = A*p;

		T alpha = rsold / ( p *Ap ); //alpha = rsold / (p' * Ap);
		
		// Update guesses
        x = x + (alpha * p);
		
		r = r - (alpha * Ap);

		// New error
		T rsnew = r*r;

		p = r + (rsnew /rsold) * p;
		
        rsold = rsnew;
		iter++;
		cout<<iter << " "<<rsold<<endl;
	}


	return x;


}


// preconditioned conjugate gradient
template <typename T>
vector<T> solveSystemIterPCG(Matrix2D<T>&A ,vector<T> &b, Matrix2D<T> & L,Matrix2D<T> & Ltrans ){
	// uses preconditioned conjugate gradient method 
	double tol = 0.1;
	// initialize guess
	vector<T> x;
	for (size_t  i = 0; i < b.size(); i++) 	x.push_back(0);

	// Find incomplete cholesky Lower matrix to precondition the solution
	//Matrix2D<T> L(b.size(),b.size());

	//this->computeIcompleteChol(L);

	// vector of residuals
	
	vector<T> r(b.size());
	
	vector<T> temp;
	//temp = this->dot(x); // A*B
	//r = b - A * x;

	r = b - (A*x);
	//for (int i = 0; i < b.size(); i++) r[i] = b[i]-temp[i]; 

	// vector Z
	vector<T> y =	forwardSubs(L,r);
	vector<T> z  = backSubs(Ltrans,y);

	// rsold = r' * r;
	T rsold=r*r;

	//for (int i = 0; i < b.size(); i++) rsold += r[i]*r[i];

	//p = r;
	vector<T> p = z;
	// iterate
	int iter =0;
	int maxIter = b.size();

	while (iter<maxIter && sqrt(rsold)>tol)
	{
		//Ap = A * p;
		//vector<double> Ap = this->dot(p);
		vector<T> Ap = A*p;
		
		//double denom=0;
		//for (int j = 0; j < b.size(); j++)  {denom += p[j] * Ap[j];}
		
		T alpha = (r*z) / (p*Ap);
		//alpha = r*z / (p' * Ap);
		//double alpha = dotVectors(r,z) / denom;

		//x = x + alpha * p;
		x = x + alpha*p;
		//for (int j = 0; j < b.size(); j++) {x[j] += alpha * p[j];}
        
		vector<T> r_old = r;
		//r = r - alpha * Ap;
		//for (int j = 0; j < b.size(); j++) {r[j] -= alpha * Ap[j];}
		r = r - alpha * Ap;
		//rsnew = r' * r;
		T rsnew =r*r;
		//double rsnew =0;
		//for (int j = 0; j < b.size(); j++) {rsnew += r[j]*r[j];}

		// z(k+1) = M\r(k+1)
		vector<T> z_old = z;
		y =	forwardSubs(L,r);
		z  = backSubs(Ltrans,r);
		//z = forwardSubs(M,r);

		T beta = (z*r) / (z_old*r_old);
		//double beta = dotVectors(z,r) / dotVectors(z_old,r_old);
		p = z + (beta * p);
		//for (int j = 0; j < b.size(); j++){ p[j] = z[j] + beta * p[j];}
		
        rsold = rsnew;
		iter++;
		cout<<iter << " "<<rsold<<endl;
	}


	return x;


}

