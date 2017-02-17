#include "ThermalFluidMesh.h"
#include "Matrix2D.h"
#include <ctime>

#define INCLUDE_MASS_MATRIX true
#define DONT_INCLUDE_MASS_MATRIX false
// Finite Elements for fluid dynamics
// Computes the sagging in an extruded pipe during cooling
// discretized with bilinear quadrilaterals.
// 


int main() {
	// double E = 210000000; // Young's modulus
	// double nu = 0.3;	  // Poisson ratio
	
	//ThermalFluidMesh FE2DMesh;
	
	//FE2DMesh.initializeFromFile("triangularMeshFine.txt");
	//FE2DMesh.initializeFromFile("triangularMesh.txt");
	//FE2DMesh.initializeFromFile("quadMesh.txt");
	//FE2DMesh.initializeFromFile("mixedMesh.txt");

	//FE2DMesh.computeElementsStifMatrix(E,nu);
	//FE2DMesh.assembleGlobalStifMatrix();
	//FE2DMesh.printDisplacements();
	
	
	ThermalFluidMesh FE2DMesh;
	FE2DMesh.initializeFromFile("pipe2Out.txt");
	clock_t timer_start =clock();
	string FN = "results/transResult10.pvd";
	ofstream ss;
	ss.open(FN);
	FE2DMesh.printPVDfileHeading(ss);
	int fileCount = 0;
	for (int step = 0; step < 1000; step++)
	{
		// print results
			if ((step % 1) == 0)
		{
			char filename[32];
			sprintf_s(filename, "results/result10_%05d.vtu", fileCount);
			FE2DMesh.printVTUfile(filename);
			FE2DMesh.prindPVDDataset(ss,fileCount,filename);
			fileCount++;
		}

		// solve for fluid and advance mesh
		FE2DMesh.assembleGlobalStifMatrix(FLUID,DONT_INCLUDE_MASS_MATRIX);
		FE2DMesh.findDisplacements(FLUID);
		FE2DMesh.advanceMesh();

		// Find temperature matrices and advance temperature field
		FE2DMesh.assembleGlobalStifMatrix(THERMAL,INCLUDE_MASS_MATRIX);
		FE2DMesh.advanceTemperatureField();
	}
	// reduce K_t
	FE2DMesh.printPVDfileFooter(ss);
	




	//FE2DMesh.assembleGlobalStifMatrix(FLUID,DONT_INCLUDE_MASS_MATRIX);
	//FE2DMesh.findDisplacements(FLUID);
	//FE2DMesh.assembleGlobalStifMatrix(FLUID,DONT_INCLUDE_MASS_MATRIX);
	//FE2DMesh.assembleTemporalMatrices();

	//FE2DMesh.findDisplacements(THERMAL);

	//FE2DMesh.assembleGlobalStifandMassMatrix(THERMAL);
	//FE2DMesh.assembleTemporalMatrices(THERMAL);
	//FE2DMesh.findDisplacements(THERMAL);

	clock_t timer_end =clock();
	FE2DMesh.printVTUfile("results9.vtu");
	
	double elapsed_secs = double(timer_end - timer_start) / CLOCKS_PER_SEC;
	
	cout<<"time elapsed " << elapsed_secs<<" s";
	
	/*
	dMatrix2D<double> A(5,5);
	A.data[0][0]= 2.;A.data[0][1]= 1.;A.data[0][2]= 1.;A.data[0][3]= 3.;A.data[0][4]= 2.; 
	A.data[1][0]= 1.;A.data[1][1]= 2.;A.data[1][2]= 2.;A.data[1][3]= 1.;A.data[1][4]= 1.; 
	A.data[2][0]= 1.;A.data[2][1]= 2.;A.data[2][2]= 9.;A.data[2][3]= 1.;A.data[2][4]= 5.; 
	A.data[3][0]= 3.;A.data[3][1]= 1.;A.data[3][2]= 1.;A.data[3][3]= 7.;A.data[3][4]= 1.; 
	A.data[4][0]= 2.;A.data[4][1]= 1.;A.data[4][2]= 5.;A.data[4][3]= 1.;A.data[4][4]= 8.; 

	vector<double> b;
	b.push_back(-2.);
	b.push_back(4.);
	b.push_back(3.);
	b.push_back(-5.);
	b.push_back(1.);

	/*
	vector<double> x = A.solveSystemLUGen(b);
	vector<double> x2 = A.solveSystemIterPCG(b);
	*/
	
	/*
	Sparserow<double> row1(10.4,5);

	row1.insertNextItem(3,7);
	row1.insertNextItem(0,6);
	row1.insertNextItem(3,2);

	row1.removeZeros();
	row1[6]=2;

	row<double> &row2 = row1;

	
	double val = row1[3];
	
	vector<double> a;
	for (int i = 0; i < 9; i++) a.push_back(2);
	double g = row1*a;


	dMatrix2D<double> A(5,5);
	row<double> & r1= A[0];
	double elValue = r1[1];

	
	A[0][0]= 2.;A[0][1]= 1.;A[0][2]= 1.;A[0][3]= 3.;A[0][4]= 2.; 
	A[1][0]= 1.;A[1][1]= 2.;A[1][2]= 2.;A[1][3]= 1.;A[1][4]= 1.; 
	A[2][0]= 1.;A[2][1]= 2.;A[2][2]= 9.;A[2][3]= 1.;A[2][4]= 5.; 
	A[3][0]= 3.;A[3][1]= 1.;A[3][2]= 1.;A[3][3]= 7.;A[3][4]= 1.; 
	A[4][0]= 2.;A[4][1]= 1.;A[4][2]= 5.;A[4][3]= 1.;A[4][4]= 8.; 
	
	spMatrix2D<double> sA(5);

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			sA.addItem(i,j,A[i][j]);
		}
	}

	vector<double> b;
	b.push_back(-2.);
	b.push_back(4.);
	b.push_back(3.);
	b.push_back(-5.);
	b.push_back(1.);

	Matrix2D<double>  & genA = A; 

	vector<double> rB = genA * b;

	vector<double> rC = matMul(A,b); //A *b;
	vector<double> rD = sA * b;
	spMatrix2D<double> spA(10);
	//row1*=2;
	
	spA.addItem(0,3,12); spA.addItem(0,4,11);  spA.addItem(0,2,13); 
	spA.addItem(1,7,12); spA.addItem(1,5,11);  spA.addItem(1,2,13);

	spA.addItem(4,4,5.3);
	spA.addItem(4,4,8.3);

	double hh = spA[1][6];
	spA[4][4]+=1000;

	spA[4][4] = A[0][0];

	spA[4][4] = 23;


	
	spA[0].removeColumn(2);

	Matrix2D<double> & A_trans = A.transpose();
	Matrix2D<double> & spA_trans = spA.transpose();

	dMatrix2D<double> M = A;
	computeIcompleteChol(M); // this has to be explicit
	M.removeZeros();
	Matrix2D<double> & Mtrans = M.transpose();
	//vector<double> backS = backSubs(sM.transpose(),b);
	//vector<double> forwardS = forwardSubs(sM,b);

	//vector<double> x =solveSystemIter( A,b);
	//vector<double> x =solveSystemIterCG( sA,b);
	//vector<double> x2 =solveSystemIterPCG( A,b,M, Mtrans);
	
	spMatrix2D<double> sC = sA-sA;
	spMatrix2D<double> sD = sA*4.0;
	spMatrix2D<double> sE = 4.0*sA;
	*/
	//cin.get();
	return 0;

}