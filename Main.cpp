#include "ThermalFluidMesh.h"
#include "Matrix2D.h"
#include <ctime>

#define INCLUDE_MASS_MATRIX true
#define DONT_INCLUDE_MASS_MATRIX false
// Finite Elements for fluid dynamics
// Computes the sagging in an extruded pipe during cooling
// discretized with biquadratic quadrilaterals.
// 


int main() {
	
	
	ThermalFluidMesh FE2DMesh;

	FE2DMesh.initializeFromFile("pipe2Out.txt");
	clock_t timer_start =clock();
	string FN = "results/transResult10.pvd";
	ofstream ss;
	ss.open(FN);
	FE2DMesh.printPVDfileHeading(ss);
	int fileCount = 0;
	for (int step = 0; step < 1000; step++)

	FE2DMesh.initializeFromFile("pipe1Out.txt");

	clock_t timer_start =clock();

	FE2DMesh.assembleGlobalStifMatrix(FLUID,false);
	FE2DMesh.findDisplacements(FLUID);
	//FE2DMesh.assembleGlobalStifandMassMatrix(THERMAL);
	//FE2DMesh.assembleTemporalMatrices(THERMAL);
	//FE2DMesh.findDisplacements(THERMAL);

	clock_t timer_end =clock();
	FE2DMesh.printVTUfile("pipeStatic.vtu");
	
	double elapsed_secs = double(timer_end - timer_start) / CLOCKS_PER_SEC;
	
	cout<<"time elapsed " << elapsed_secs<<" s";
	
	/*
	dMatrix2D<double> A(5,5);
	A.data[0][0]= 2.;A.data[0][1]= 1.;A.data[0][2]= 1.;A.data[0][3]= 3.;A.data[0][4]= 2.; 
	A.data[1][0]= 1.;A.data[1][1]= 2.;A.data[1][2]= 2.;A.data[1][3]= 1.;A.data[1][4]= 1.; 
	A.data[2][0]= 1.;A.data[2][1]= 2.;A.data[2][2]= 9.;A.data[2][3]= 1.;A.data[2][4]= 5.; 
	A.data[3][0]= 3.;A.data[3][1]= 1.;A.data[3][2]= 1.;A.data[3][3]= 7.;A.data[3][4]= 1.; 
	A.data[4][0]= 2.;A.data[4][1]= 1.;A.data[4][2]= 5.;A.data[4][3]= 1.;A.data[4][4]= 8.; 

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

	clock_t timer_end =clock();

	FE2DMesh.printVTUfile("results9.vtu");
	
	double elapsed_secs = double(timer_end - timer_start) / CLOCKS_PER_SEC;
	
	cout<<"time elapsed " << elapsed_secs<<" s";
	
	return 0;

}