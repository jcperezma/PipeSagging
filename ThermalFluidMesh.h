#pragma once
#include "MathUtils.h"
#include "FiniteElement.h"
#include <memory>
#include <sstream> 
#include <fstream>
#include <string>
#include <iostream>



class ThermalFluidMesh {

public:
	virtual void initializeFromFile(const string & fileName);
	virtual void assembleGlobalStifMatrix(physType const & type, bool computeM);
	virtual void assembleTemporalMatrices(physType const & type);
	virtual void findDisplacements(physType const & type);
	virtual void advanceMesh();

	// output operations
	virtual void prindPVDDataset (ofstream &out, int step, string filename );
	virtual void printPVDfileHeading(ofstream &out);
	virtual void printPVDfileFooter(ofstream & out);
	virtual void printVTUfile(string fileName);
	virtual void printDisplacements();

protected:
	vector<unique_ptr<FiniteElement>> elements;
	vector<Point2D> coords;
	SimParameters params;   // Values needed in the simulation, material properties etc...
	// Arrays for Velocities
	FE_data fluidData;

	// Arrays for Temperature
	FE_data thermalData;


};