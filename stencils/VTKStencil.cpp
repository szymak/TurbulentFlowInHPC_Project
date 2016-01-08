#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "VTKStencil.h"
#include "Iterators.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>

VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
	this->_outputFile = new std::ofstream;
	this->_velocityStringStream << std::fixed << std::setprecision(6);
	this->_pressureStringStream << std::fixed << std::setprecision(6);
}

VTKStencil::~VTKStencil () {
	delete this->_outputFile;
}

void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
	
	// Skip the ghost cell
	if(i == 1 || j == 1) return;
	
    const int obstacle = flowField.getFlags().getValue(i, j);
    
	if(obstacle & OBSTACLE_SELF) {
		this->_pressureStringStream << "0.0" << std::endl;
		this->_velocityStringStream << "0.0 0.0 0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[2];
		flowField.getPressureAndVelocity(pressure, velocity, i, j);
		this->_pressureStringStream << pressure << std::endl;
		this->_velocityStringStream << velocity[0] << " " << velocity[1] << " 0.0" << std::endl;
	}
	
}

void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {

	// Skip the ghost cell
	if(i == 1 || j == 1 || k == 1) return;
	
    const int obstacle = flowField.getFlags().getValue(i, j, k);
    
	if(obstacle & OBSTACLE_SELF) {
		this->_pressureStringStream << "0.0" << std::endl;
		this->_velocityStringStream << "0.0 0.0 0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[3];
		flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
		this->_pressureStringStream << pressure << std::endl;
		this->_velocityStringStream << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
	}
	
}

void VTKStencil::write ( FlowField & flowField, int timeStep, std::string foldername ) {

	std::cout << "=== Writing VTK Output ===" << std::endl;
	
	// Open the file and set precision
	this->_outputFile->open(this->getFilename(timeStep, foldername).c_str());
	*this->_outputFile << std::fixed << std::setprecision(6);
	
	// Output the different sections of the file
	this->writeFileHeader();
	this->writeGrid(flowField);
	this->writeCellDataHeader(flowField);
	this->writePressure();
	this->writeVelocity();
	
	// Close the file
	this->_outputFile->close();
	
	this->clearStringStreams();
	
}

void VTKStencil::writeFileHeader() {
	*this->_outputFile << "# vtk DataFile Version 2.0" << std::endl;
	*this->_outputFile << "Turbulent Flow Simulation" << std::endl;
	*this->_outputFile << "ASCII" << std::endl << std::endl;
}

void VTKStencil::writeGrid ( FlowField & flowField ) {

	int is3D = this->_parameters.geometry.dim == 3;
	
	// Get number of grid points
	int gridX = flowField.getNx() + 1;
	int gridY = flowField.getNy() + 1;
	int gridZ = is3D ? flowField.getNz() + 1 : 1;
	
	// Print grid header
	*this->_outputFile << "DATASET STRUCTURED_GRID" << std::endl;
	*this->_outputFile << "DIMENSIONS " << gridX << " " << gridY  << " " << gridZ << std::endl;
	*this->_outputFile << "POINTS " << gridX * gridY * gridZ << " float" << std::endl;
	
	// Get id of first and last cell in each dimension
	int firstCell = 2; // 0 and 1 are ghost cells
	int lastCellX = flowField.getCellsX();
	int lastCellY = flowField.getCellsY();
	int lastCellZ = is3D ? flowField.getCellsZ() : 3;
	
	// Print grid
	for(int k = firstCell; k < lastCellZ; k++) {
		for(int j = firstCell; j < lastCellY; j++) {
			for(int i = firstCell; i < lastCellX; i++) {
				*this->_outputFile << this->_parameters.meshsize->getPosX(i, j, k) << " " << this->_parameters.meshsize->getPosY(i, j, k) << " " << this->_parameters.meshsize->getPosZ(i, j, k) << std::endl;
			}
		}
	}
	
	*this->_outputFile << std::endl;
	
}

void VTKStencil::writeCellDataHeader ( FlowField & flowField ) {

	int is3D = this->_parameters.geometry.dim == 3;
	int cellsX = flowField.getNx();
	int cellsY = flowField.getNy();
	int cellsZ = flowField.getNz();
	int numCells = is3D ? (cellsX * cellsY * cellsZ) : (cellsX * cellsY);
	*this->_outputFile << "CELL_DATA " << numCells << std::endl;
	
}

void VTKStencil::writePressure ( ) {

	// Print header
	*this->_outputFile << "SCALARS pressure float 1" << std::endl;
	*this->_outputFile << "LOOKUP_TABLE default" << std::endl;
	
	// Print pressure values
	*this->_outputFile << this->_pressureStringStream.str().c_str();
	
	*this->_outputFile << std::endl;
	
}

void VTKStencil::writeVelocity ( ) {

	// Print header
	*this->_outputFile << "VECTORS velocity float" << std::endl;
	
	// Print velocity values
	*this->_outputFile << this->_velocityStringStream.str().c_str();
	
	*this->_outputFile << std::endl;
	
}

std::string VTKStencil::getFilename( int timeStep, std::string foldername ) {	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::stringstream filename;
	filename << foldername << "/" << this->_parameters.vtk.prefix << "." << rank << "." << timeStep << ".vtk";
	return filename.str();
}

void VTKStencil::clearStringStreams() {
	_pressureStringStream.str("");
	_velocityStringStream.str("");
}
