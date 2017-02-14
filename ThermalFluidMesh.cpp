#include "ThermalFluidMesh.h"



void ThermalFluidMesh::assembleGlobalStifMatrix(physType const & type, bool computeM){

	// compute stiffness Matrix for each element and add it to the global stiffness matrix
	// reset matrix K and vector F;
	
	dMatrix2D<double> k;
	vector<double> f;
	dMatrix2D<double> m;
	switch (type)
	{
	case FLUID:
		k.resize(18,18);
		f.resize(18);
		fluidData.K.setToZero();
		for (int i = 0; i < fluidData.F.size(); i++) fluidData.F[i] =0;
		if (computeM){
			m.resize(18,18);
			fluidData.M.setToZero();
		}

	break;
	case THERMAL:
		k.resize(9,9);
		f.resize(9);
		thermalData.K.setToZero();
		for (int i = 0; i < thermalData.F.size(); i++) thermalData.F[i] =0;
		if (computeM){
			m.resize(9,9);
			thermalData.M.setToZero();
		}
	break;
	default:
		break;
	}
	
	for (auto &element: elements )
	{
		element->computeLocalStifandMassMatrix(coords,thermalData.U,fluidData.U,k,m,f,type,params,computeM);

		switch (type)
		{
		case FLUID:
			element->addLocalStifandMassToGlobal(k,m,f,fluidData,type,computeM);
				break;
		case THERMAL:
			element->addLocalStifandMassToGlobal(k,m,f,thermalData,type,computeM);
				break;
		default:
			break;
		}
	}
}

void ThermalFluidMesh::advanceMesh(){
	// advance the coords based on the velocities

	for (int i = 0; i < coords.size(); i++)
	{
		coords[i].x += params.dt*fluidData.U[i*fluidData.dim];
		coords[i].y += params.dt*fluidData.U[i*fluidData.dim+1];

	}

};



void printMatrix(string fileName, dMatrix2D<double> mat){
	ofstream out(fileName);
	for (int i = 0; i < mat.getnRows(); i++)
	{
		for (int j = 0; j < mat.getNColumns(); j++)
		{
			out<< mat[i][j] << " ";
		}
		out<<endl;
	}
	out.close();
}


void printVector(string fileName, vector<double> vec){

	ofstream out(fileName);
	for (int i = 0; i < vec.size(); i++)
	{

			out<< vec[i]<< endl;

	}
	out.close();

}

void ThermalFluidMesh::prindPVDDataset (ofstream &out, int step, string filename ){
	out << "<DataSet timestep=\"" << step << "\"  group = \"\" part=\"0\" file= \""<< filename << "\"/>" <<endl;
	

}

void ThermalFluidMesh::printPVDfileHeading(ofstream &out){
	
	// header
	out << "<VTKFile type=\"" << "Collection" << "\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	
	out << "  <Collection>" << endl;
	

}

void ThermalFluidMesh::printPVDfileFooter(ofstream & out){
	
	// footer
	out << "  </Collection>" << endl;
	out << "</VTKFile>";

}

void ThermalFluidMesh::assembleTemporalMatrices(physType const & type){
	// The idea is to end up with a system of equations
	// of the type K_t T_new = F_t
	// where 
	// F_t = ( M / dt -(1-theta) * K ) T_old
	// K_t = ( M / dt + theta * K )
	
	vector<double> F_t;
	spMatrix2D<double>K_t = thermalData.M * (1/params.dt ) + params.theta * thermalData.K;
	//initialize K_reduced, depends if it is dense or sparse
	// this is actually the one that has to be sparse
	spMatrix2D<double> K_treduced;
	K_t.reduce(K_treduced,thermalData.DOFid,thermalData.BC_U_id); // did it this way to have the definition of K and K_reduced in the same place

	// apply boundary conditions to F_t and reduce F_t

	// 1) add body force boundary conditions to right hand side

	for (int i = 0; i < thermalData.BC_F.size(); i++) thermalData.F[thermalData.BC_F_id[i]] = thermalData.BC_F[i];

	string FN = "transResult.pvd";
	ofstream ss;
	ss.open(FN);
	
	printPVDfileHeading(ss);

	int fileCount = 0;
	for (int step = 0; step < 500; step++)
	{
		if ((step % 20) == 0)
	{
		char filename[24];
		sprintf_s(filename, "result_%04d.vtu", fileCount);
		printVTUfile(filename);
		prindPVDDataset(ss,step,filename);
		fileCount++;
	}

		F_t= (thermalData.M * (1/params.dt) - (1-params.theta)* thermalData.K ) * thermalData.U; 

	// 2)Add displacements boundary conditions
	//  modify right hand side  vector to account for the velocity BC
	int num_displacement_DOF = thermalData.DOFid.size();
	for (int i = 0; i < thermalData.BC_U_id.size();i++)
	{
		for (int j = 0; j < num_displacement_DOF; j++)
		{
			int rowIndex = thermalData.DOFid[j];
			int columnIndex = thermalData.BC_U_id[i];
			F_t[rowIndex] -= K_t[rowIndex][columnIndex]*thermalData.BC_U[i];
		}
	}

	// 3)create reduced right hand side vector 
	vector<double>F_reduced(num_displacement_DOF); 
	for (int i = 0; i < num_displacement_DOF; i++)
	{
		F_reduced[i] = F_t[thermalData.DOFid[i]];
	}

	// 4)create reduced stiffnes matrix
	// remove columns and rows of the global stiffness matrix of the nodal displacements that were given by the BC.

	

	//printMatrix("matrix.txt",data.K_reduced);
	

	// 5)solve the reduced system of equations
	//vector<double> reduced_u = solveSystemIter(K_reduced,F_reduced);
	//vector<double> reduced_u2 = solveSystemIterOld(K_reduced,F_reduced);
	vector<double> reduced_u = solveSystemIterCG(K_treduced,F_reduced,thermalData.u_reduced); //U_new
	//vector<double> reduced_u = K_reduced.solveSystemIterPCG(F_reduced);


	// 6)add the reduced displacements to the global array of displacements.
	for (int i = 0; i < num_displacement_DOF; i++)
	{
		thermalData.U[thermalData.DOFid[i]] = reduced_u[i];
	}
		
	
			
	}
	// reduce K_t
	printPVDfileFooter(ss);
	// solve for T_new 

	
}

void DoFindDisplacements( FE_data & data  ){
	
	// Generalized way to solve a finite element system once the Global stiffness matrix K and Loads vector F have been assembled

	
	// 1) add body force boundary conditions to right hand side

	for (int i = 0; i < data.BC_F.size(); i++) data.F[data.BC_F_id[i]] = data.BC_F[i];

	
	// 2)Add displacements boundary conditions
	//  modify right hand side  vector to account for the velocity BC
	int num_displacement_DOF = data.DOFid.size();
	for (int i = 0; i < data.BC_U_id.size();i++)
	{
		for (int j = 0; j < num_displacement_DOF; j++)
		{
			int rowIndex = data.DOFid[j];
			int columnIndex = data.BC_U_id[i];
			data.F[rowIndex] -= data.K[rowIndex][columnIndex]*data.BC_U[i];
		}
	}

	// 3)create reduced right hand side vector 
	vector<double>F_reduced(num_displacement_DOF); 
	for (int i = 0; i < num_displacement_DOF; i++)
	{
		F_reduced[i] = data.F[data.DOFid[i]];
	}

	// 4)create reduced stiffnes matrix
	// remove columns and rows of the global stiffness matrix of the nodal displacements that were given by the BC.

	//initialize K_reduced, depends if it is dense or sparse
	// this is actually the one that has to be sparse

	data.K.reduce(data.K_reduced,data.DOFid,data.BC_U_id); // did it this way to have the definition of K and K_reduced in the same place


	//printMatrix("matrix.txt",data.K_reduced);
	printVector("vector.txt",F_reduced);

	// 5)solve the reduced system of equations
	//vector<double> reduced_u = solveSystemIter(K_reduced,F_reduced);
	//vector<double> reduced_u2 = solveSystemIterOld(K_reduced,F_reduced);
	data.u_reduced = solveSystemIterCG(data.K_reduced,F_reduced,data.u_reduced);
	//vector<double> reduced_u = K_reduced.solveSystemIterPCG(F_reduced);
	//vector<double> reduced_u = K_reduced.solveSystemLUGen(F_reduced);
	//vector<double> reduced_u = K_reduced.solveSystemLUChol(F_reduced);
	printVector("u.txt",data.u_reduced);

	// 6)add the reduced displacements to the global array of displacements.
	for (int i = 0; i < num_displacement_DOF; i++)
	{
		data.U[data.DOFid[i]] = data.u_reduced[i];
	}

}


void ThermalFluidMesh::findDisplacements(physType const & type){

	switch (type)
	{
	case FLUID:
		DoFindDisplacements(fluidData);

		break;
	case THERMAL:
		DoFindDisplacements(thermalData);
		break;

	default:
		break;
	}
	
	
}

void writeVtkVectorArray(ostream &out, const string &name, const vector<Point2D> &data) {
	out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;
	for (int c = 0; c < data.size(); c++) {
		
		out <<"			 "<< data[c].x << " " << data[c].y <<" " << 0 <<endl;
	}
	out << "		 </DataArray>" << endl;
}

void writeVtkVectorArrayVelocity(ostream &out,const string &name,vector<double>const  & U, int numNodes) {
	out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" Format=\"ascii\">"<<endl;
	for (int c = 0; c < numNodes; c++) {
		
		out <<"			 "<<  U[c*2] << " " << U[c*2+1] <<" " <<0<<endl;
	}
	out << "	   </DataArray>" << endl;
}

void writeVtkVectorArrayTemp(ostream &out,const string &name,vector<double>const  & T) {
	out << "        <DataArray type=\"Float64\" Name=\"" << name << " \" format=\"ascii\">";
	for (int c = 0; c < T.size(); c++) {
		
		out << " " << T[c] ;
	}
	out << "    </DataArray>" << endl;
}

void writeVtkCellData(ostream &out, int numPoints, double temperature) {
	out << "      <CellData Scalars=\"Temp\">" << endl;
	out << "        <DataArray type=\"Float64\" Name=\"Temp\" format=\"ascii\">";
	for (int c = 0; c < numPoints; c++) {
		out << " " << temperature;
	}
	out << " </DataArray>" << endl;
	out << "      </CellData>" << endl;
}



void ThermalFluidMesh::printVTUfile(string fileName){
	ofstream out(fileName);
	// header
	out << "<VTKFile type=\"" << "UnstructuredGrid" << "\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	
	out << "  <UnstructuredGrid>" << endl;
	out << "    <Piece NumberOfPoints=\"" << coords.size() << "\" NumberOfCells=\"" << elements.size() << "\" >" << endl;
	
	out << "      <Points>" << endl;
	writeVtkVectorArray(out, "Positions", coords);
	out << "      </Points>" << endl;
	
	out << "      <Cells>" << endl; {
		out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<endl; {
			for (auto &elem : this->elements) {
				out<<"			";
				for (auto &index : elem->getNodes()) {
					out << " " << index;
				}
				out<<endl;
			}
		} out << "		</DataArray>" << endl;

		out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"; {
			int endIndex = 0;
			for (auto &elem : this->elements) {
				endIndex += 9;
				out << " " << endIndex;
			}
		} out << " </DataArray>" << endl;

		out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"; {
			for (auto &elem : this->elements) {
				out << " " << 23;
			}
		} out << " </DataArray>" << endl;
	} out << "      </Cells>" << endl;

	// write values
	out << "      <PointData Scalars= \"scalars\" \>" << endl;
	if (fluidData.U.size()>0) writeVtkVectorArrayVelocity(out,"Velocity",fluidData.U,coords.size());
	if (thermalData.U.size()>0) writeVtkVectorArrayTemp(out,"Temperature",thermalData.U);
	out << "      </PointData>" << endl;

	//writeVtkCellData(out,coords.size(),45);
	//writeVtkVectorArray(out, "Velocity", time.velocities);
	//writeVtkVectorArray(out, "Vorticity", time.vorticities);
	

	out << "    </Piece>" << endl;
	out << "  </UnstructuredGrid>" << endl;
	out << "</VTKFile>" << endl;

}

void readNodesCoords(ifstream & tempFile,vector<Point2D>&coords){
	stringstream ss;
	string line;
	// read number of nodes

	getline (tempFile,line);
	ss<<line;
	int NumNodes;
	ss>>NumNodes;
	ss=stringstream();// flush stream

	//read node coords
	for (int i = 0; i < NumNodes; i++)
	{
		Point2D p;
		getline (tempFile,line);
		ss<<line;
		ss>>p.x;
		ss>>p.y;
		ss=stringstream();// flush the stream
		coords.push_back(p);
	}
}

void readElements(ifstream & tempFile,vector<unique_ptr<FiniteElement>>&elements){
	stringstream ss;
	string line;
	getline (tempFile,line);
	ss<<line;
	int NumElements;
	ss>>NumElements;
	ss=stringstream();// flush stream

	// read element topology
	for (int i = 0; i < NumElements; i++)
	{
		vector<int> nodeID;
		int node;
		getline (tempFile,line);
		ss<<line;
		
		for (int j = 0; j < 9; j++)
		{
			ss>>node;
			nodeID.push_back(node);
		}
		unique_ptr<FiniteElement> BLiQuadElement (new BiLinearQuadThermalFluid);
		BLiQuadElement->setNodeIds(nodeID);
		elements.push_back(move(BLiQuadElement));

		ss=stringstream();// flush the stream
	}
}

void readBC(ifstream & tempFile,FE_data &FE_vars){
	stringstream ss;
	string line;
	// read displacement constrains
	getline (tempFile,line);
	ss<<line;
	int NumUConstrains;
	ss>>NumUConstrains;
	ss=stringstream();// flush stream

	//read indices of the constrains
	for (int i = 0; i < NumUConstrains; i++)
	{
		int U_ID;
		getline (tempFile,line);
		ss<<line;
		ss>>U_ID;
		ss=stringstream();// flush stream
		FE_vars.BC_U_id.push_back(U_ID);
	}

	// read displacement constrain values
	for (int i = 0; i < NumUConstrains; i++)
	{
		double U_value;
		getline (tempFile,line);
		ss<<line;
		ss>>U_value;
		ss=stringstream();// flush stream
		FE_vars.BC_U.push_back(U_value);
	}

	// read nodal forces
	getline (tempFile,line);
	ss<<line;
	int NumNodalForces;
	ss>>NumNodalForces;
	ss=stringstream();// flush stream

	//read indices of the constrains
	for (int i = 0; i < NumNodalForces; i++)
	{
		int F_ID;
		getline (tempFile,line);
		ss<<line;
		ss>>F_ID;
		ss=stringstream();// flush stream
		FE_vars.BC_F_id.push_back(F_ID);
	}

	// read displacement constrain values
	for (int i = 0; i < NumNodalForces; i++)
	{
		double F_value;
		getline (tempFile,line);
		ss<<line;
		ss>>F_value;
		ss=stringstream();// flush stream
		FE_vars.BC_F.push_back(F_value);
	}
}


void initK(dMatrix2D<double>& K, vector<unique_ptr<FiniteElement>>&elements, int NumNodes, int dim ){
	K.resize(NumNodes*dim,NumNodes*dim);
	K =0;
}


void initK(spMatrix2D<double>& K, spMatrix2D<double>& M, vector<unique_ptr<FiniteElement>>& elements, int NumNodes, int dim ){
// For the sparse matrix case

	// First resize Matrix
	K.resize(NumNodes*dim,NumNodes);
	M.resize(NumNodes*dim,NumNodes);
	// Loop through the elements 

	for (auto const & element : elements)
	{
		for (auto const & a : element->getNodes())
		{
			int p = dim * a;
			for (auto const & b : element->getNodes())
			{
				int q = dim *b;
				for (int i = 0; i < dim; i++)
				{
					for (int j = 0; j < dim; j++)
					{
						K.addItem(p+i,q+j,0.);
						M.addItem(p+i,q+j,0.);
					}
				}
			}
		}
	}

}

void initGlobals(FE_data &Fe_vars, vector<unique_ptr<FiniteElement>>&elements, int NumNodes ){
	
	initK(Fe_vars.K,Fe_vars.M,elements,NumNodes,Fe_vars.dim);

	// allocate and initialize displacementes array
	Fe_vars.U.resize(NumNodes*Fe_vars.dim);
	for (auto displacement : Fe_vars.U)
	{
		displacement =0.0;
	}

	// add prescribed displacements to global array of displacements
	for (int i=0; i< Fe_vars.BC_U_id.size(); i++) {Fe_vars.U[Fe_vars.BC_U_id[i]] = Fe_vars.BC_U[i];}

	// allocate and initialize global force vectors
	Fe_vars.F.resize(NumNodes*Fe_vars.dim);
	for (auto force : Fe_vars.F)
	{
		force =0.0;
	}

	// initialize array with the indices to the degrees of fredom 
	int count =0;

	// for the velocity

	for (int i = 0; i < Fe_vars.U.size(); i++)
	{
		if (count < Fe_vars.BC_U_id.size() && i != Fe_vars.BC_U_id[count])
		{
			// it has not been prescribed
			Fe_vars.DOFid.push_back(i);
		}else
		{
			// it has been prescribed check next prescribed index
			count ++;
		}
	}
// initialize reduced U, its size is the number of degrees of freedom

	int numDOF = Fe_vars.DOFid.size();
	Fe_vars.u_reduced.resize(numDOF);
	for (int i = 0; i < numDOF; i++)
	{
		Fe_vars.u_reduced[i] = 0;
	}

}




void ThermalFluidMesh::initializeFromFile(const string & fileName){
	// the format of the input file is:
	// rho mu lambda g_y nGauss
	// numNodes
	// x1 y1
	//  . .
	// xN yN
	// numElements
	// indexNode1 indexNode2 ... // there's just one type of elements.
	// . .. .. . .. .
	// indexNode1 indexNode2 ..
	// numVelocityConstrians
	// displacementID1
	// ...
	// displacemnteIDN
	// velocityValue1
	// ...
	// velocityValueN
	// numVelocityLoads
	// nodalforceID1
	// ...
	// nodalforceIDN
	// nodalforceValue1
	// ...
	// nodalForceValueN
	// numTemperatureConstrians
	// TemperatureID1
	// ...
	// TemperatureIDN
	// TemperatureValue1
	// ...
	// TemperatureValueN
	// numVelocityLoads
	// tempforceID1
	// ...
	// tempforceIDN
	// tempforceValue1
	// ...
	// tempForceValueN


	string line;
	stringstream ss;
	//open file
	ifstream tempFile (fileName);
	bool didntFail = tempFile.is_open();

	// read simulation parameters
	getline (tempFile,line);
	ss<<line;
	ss>>this->params.rho; ss>>this->params.mu; ss>>this->params.lambda; ss>>this->params.g_y; ss>>this->params.nGauss;
	ss>>this->params.kappa; ss>>params.theta; ss>>params.dt; ss>>params.Cp;ss>>thermalData.U_o;
	ss=stringstream();// flush stream

	// read mesh topology
	readNodesCoords(tempFile,coords); // read number of nodes and coords
	
	// read element connectivity
	readElements(tempFile,elements);

	// read BC for the velocity
	readBC(tempFile, fluidData);
	
	// read BC for the Temperature
	
	// this depends if it is transcient
	readBC(tempFile, thermalData);

	tempFile.close();
	// initialize global stiffness matrix, and other data structures, for sparse matrix the connectivity is needed. 
	
	// For the fluid variables
	fluidData.dim=2;
	initGlobals(fluidData, elements, coords.size() );
	

	// for the thermal variables
	thermalData.dim=1;
	initGlobals(thermalData, elements, coords.size() );

	//Transcient
	//thermalData_trans.dim=1;
	//initGlobals(thermalData_trans, elements, coords.size() );
	

}


void ThermalFluidMesh::printDisplacements(){
	cout<<"Displacements:\n\n";
	cout<<"DOF_ID	Displacement [m] \n\n";
	int DOF_ID=0;
	for (auto displacement : fluidData.U)
	{
		cout<<DOF_ID <<	"	 "<< displacement<<"\n";
		DOF_ID++;
	}
}