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
	virtual void printDisplacements();
	virtual void printVTUfile(string fileName);
	virtual void advanceMesh();
protected:
	virtual void DoAssembleGlobalStifMatrix(physType const & type, FE_data & FEdata, bool computeM);

	vector<unique_ptr<FiniteElement>> elements;
	vector<Point2D> coords;
	SimParameters params;   // Values needed in the simulation, material properties etc...
	// Arrays for Velocities
	FE_data fluidData;

	// Arrays for Temperature
	FE_data thermalData;
	

};