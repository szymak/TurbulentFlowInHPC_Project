#include <iostream>
#include <sstream>
#include <iomanip>
#include "TurbulentVTKStencil.h"
#include "StencilFunctions.h"

TurbulentVTKStencil::TurbulentVTKStencil(Parameters & parameters) :
	VTKStencil(parameters),
	FieldStencil<TurbulentFlowField>(parameters)
{
}

void TurbulentVTKStencil::apply(TurbulentFlowField & flowField, int i, int j) {

	const int obstacle = flowField.getFlags().getValue(i, j);

	if(obstacle & OBSTACLE_SELF) {
		_pressureStringStream << "0.0" << std::endl;
		_velocityStringStream << "0.0 0.0 0.0" << std::endl;
		_turbulentViscosityStringStream << "0.0" << std::endl;
		_viscousStressXYStringStream << "0.0" << std::endl;
		_ReStressXYStringStream << "0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[2];
		flowField.getPressureAndVelocity(pressure, velocity, i, j);
		computeShearStresses(flowField, i, j);
		_pressureStringStream << pressure << std::endl;
		_velocityStringStream << velocity[0] << " " << velocity[1] << " 0.0" << std::endl;
		_turbulentViscosityStringStream << flowField.getTurbulentViscosity().getScalar(i, j) << std::endl;
		_viscousStressXYStringStream << _shearStresses[0] / VTKStencil::_parameters.flow.Re << std::endl;
		_ReStressXYStringStream << _shearStresses[0] * flowField.getTurbulentViscosity().getScalar(i, j) << std::endl;
	}
	
}

void TurbulentVTKStencil::apply(TurbulentFlowField & flowField, int i, int j, int k) {

	const int obstacle = flowField.getFlags().getValue(i, j, k);

	if(obstacle & OBSTACLE_SELF) {
		_pressureStringStream << "0.0" << std::endl;
		_velocityStringStream << "0.0 0.0 0.0" << std::endl;
		_turbulentViscosityStringStream << "0.0" << std::endl;
		_viscousStressXYStringStream << "0.0" << std::endl;
		_viscousStressXZStringStream << "0.0" << std::endl;
		_ReStressXYStringStream << "0.0" << std::endl;
		_ReStressXZStringStream << "0.0" << std::endl;
	} else {
		FLOAT pressure;
		FLOAT velocity[3];
		flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
		computeShearStresses(flowField, i, j, k);
		_pressureStringStream << pressure << std::endl;
		_velocityStringStream << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
		_turbulentViscosityStringStream << flowField.getTurbulentViscosity().getScalar(i, j, k) << std::endl;
		_viscousStressXYStringStream << _shearStresses[0] / VTKStencil::_parameters.flow.Re << std::endl;
		_viscousStressXZStringStream << _shearStresses[1] / VTKStencil::_parameters.flow.Re << std::endl;
		_ReStressXYStringStream << _shearStresses[0] * flowField.getTurbulentViscosity().getScalar(i, j, k) << std::endl;
		_ReStressXZStringStream << _shearStresses[1] * flowField.getTurbulentViscosity().getScalar(i, j, k) << std::endl;
	}
	
}


void TurbulentVTKStencil::write(TurbulentFlowField & flowField, int timeStep, std::string foldername) {

	std::cout << "=== Writing VTK Output ===" << std::endl;

	// Open the file and set precision
	_outputFile->open(getFilename(timeStep, foldername).c_str());
	*_outputFile << std::fixed << std::setprecision(6);

	// Output the different sections of the file
	writeFileHeader();
	writeGrid(flowField);
	writeCellDataHeader(flowField);
	writePressure();
	writeVelocity();
	writeTurbulentViscosity();
	writeShearStresses();

	// Close the file
	_outputFile->close();

	// Clear string streams
	clearStringStreams();
	
}

void TurbulentVTKStencil::writeTurbulentViscosity() {

	// Print header
	*_outputFile << "SCALARS turbulent_viscosity float 1" << std::endl;
	*_outputFile << "LOOKUP_TABLE default" << std::endl;

	// Print turbulent viscosity values
	*_outputFile << _turbulentViscosityStringStream.str().c_str();

	*_outputFile << std::endl;

}

void TurbulentVTKStencil::writeShearStresses() {

	*_outputFile << "SCALARS viscous_stress_XY float 1" << std::endl;
	*_outputFile << "LOOKUP_TABLE default" << std::endl;
	*_outputFile << _viscousStressXYStringStream.str().c_str();

	*_outputFile << "SCALARS Re_stress_XY float 1" << std::endl;
	*_outputFile << "LOOKUP_TABLE default" << std::endl;
	*_outputFile << _ReStressXYStringStream.str().c_str();
	
	if(VTKStencil::_parameters.geometry.dim == 3) {
		*_outputFile << "SCALARS viscous_stress_XZ float 1" << std::endl;
		*_outputFile << "LOOKUP_TABLE default" << std::endl;
		*_outputFile << _viscousStressXZStringStream.str().c_str();
	
		*_outputFile << "SCALARS Re_stress_XZ float 1" << std::endl;
		*_outputFile << "LOOKUP_TABLE default" << std::endl;
		*_outputFile << _ReStressXZStringStream.str().c_str();
	}
	
	*_outputFile << std::endl;
	
}


void TurbulentVTKStencil::computeShearStresses(TurbulentFlowField & flowField, int i, int j) {
	loadLocalVelocity2D(flowField, _localVelocity, i, j);
	loadLocalMeshsize2D(VTKStencil::_parameters, _localMeshsize, i, j);
	_shearStresses[0] = dudy(_localVelocity, _localMeshsize);
}

void TurbulentVTKStencil::computeShearStresses(TurbulentFlowField & flowField, int i, int j, int k) {
	loadLocalVelocity3D(flowField, _localVelocity, i, j, k);
	loadLocalMeshsize3D(VTKStencil::_parameters, _localMeshsize, i, j, k);
	_shearStresses[0] = dudy(_localVelocity, _localMeshsize);
	_shearStresses[1] = dudz(_localVelocity, _localMeshsize);
}


void TurbulentVTKStencil::clearStringStreams() {
	_pressureStringStream.str("");
	_velocityStringStream.str("");
	_turbulentViscosityStringStream.str("");
	_viscousStressXYStringStream.str("");
	_viscousStressXZStringStream.str("");
	_ReStressXYStringStream.str("");
	_ReStressXZStringStream.str("");
}


