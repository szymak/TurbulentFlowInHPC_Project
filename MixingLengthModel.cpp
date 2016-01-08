#include "MixingLengthModel.h"
#include <iostream>


MixingLengthModel::MixingLengthModel(const Parameters & parameters) :
	_parameters(parameters)
{	
}
MixingLengthModel::~MixingLengthModel() {}

/**************************************************************************************
	Lm = Kh 
**************************************************************************************/

LmKh::LmKh(const Parameters & parameters) :
	MixingLengthModel(parameters)
{
}

FLOAT LmKh::at(TurbulentFlowField & flowField, int i, int j) {
	return 0.41 * flowField.getWallDistance(i, j);
}

FLOAT LmKh::at(TurbulentFlowField & flowField, int i, int j, int k) {
	return 0.41 * flowField.getWallDistance(i, j, k);
}


/**************************************************************************************
	Lm = min(Kh, 0.09 delta) 
	Assuming laminar flat plate Blasius boundary layer
**************************************************************************************/

LmLaminarFlatPlate::LmLaminarFlatPlate(const Parameters & parameters) :
	MixingLengthModel(parameters)
{
}

FLOAT LmLaminarFlatPlate::at(TurbulentFlowField & flowField, int i, int j) {
	FLOAT kh = 0.41 * flowField.getWallDistance(i, j);
	FLOAT x = _parameters.meshsize->getPosX(i, j) + 0.5 * _parameters.meshsize->getDx(i, j);
	FLOAT centerLineVelocity = flowField.getCenterLineVelocity()[i];
	// FLOAT centerLineVelocity = flowField.getVelocity().getVector(i, _parameters.geometry.sizeY / 2)[0];
	FLOAT downstreamRe = centerLineVelocity * x * _parameters.flow.Re;
	FLOAT delta = 4.91 * x / pow(downstreamRe, 0.5);
	return std::min(kh, 0.09 * delta);
}

FLOAT LmLaminarFlatPlate::at(TurbulentFlowField & flowField, int i, int j, int k) {
	FLOAT kh = 0.41 * flowField.getWallDistance(i, j, k);
	FLOAT x = _parameters.meshsize->getPosX(i, j, k) + 0.5 * _parameters.meshsize->getDx(i, j, k);
	FLOAT centerLineVelocity = flowField.getCenterLineVelocity()[i];
	// FLOAT centerLineVelocity = flowField.getVelocity().getVector(i, _parameters.geometry.sizeY / 2, _parameters.geometry.sizeZ / 2)[0];
	FLOAT downstreamRe = centerLineVelocity * x * _parameters.flow.Re;
	FLOAT delta = 4.91 * x / pow(downstreamRe, 0.5);
	return std::min(kh, 0.09 * delta);
}


/**************************************************************************************
	Lm = min(Kh, 0.09 delta) 
	Assuming turbulent flat plate Blasius boundary layer
**************************************************************************************/

LmTurbulentFlatPlate::LmTurbulentFlatPlate(const Parameters & parameters) :
	MixingLengthModel(parameters)
{
}

FLOAT LmTurbulentFlatPlate::at(TurbulentFlowField & flowField, int i, int j) {
	FLOAT kh = 0.41 * flowField.getWallDistance(i, j);
	FLOAT x = _parameters.meshsize->getPosX(i, j) + 0.5 * _parameters.meshsize->getDx(i, j);
	FLOAT centerLineVelocity = flowField.getCenterLineVelocity()[i];
	// FLOAT centerLineVelocity = flowField.getVelocity().getVector(i, _parameters.geometry.sizeY / 2)[0];
	FLOAT downstreamRe = centerLineVelocity * x * _parameters.flow.Re;
	FLOAT delta = 0.382 * x / pow(downstreamRe, 0.2);
	return std::min(kh, 0.09 * delta);
}

FLOAT LmTurbulentFlatPlate::at(TurbulentFlowField & flowField, int i, int j, int k) {
	FLOAT kh = 0.41 * flowField.getWallDistance(i, j, k);
	FLOAT x = _parameters.meshsize->getPosX(i, j, k) + 0.5 * _parameters.meshsize->getDx(i, j, k);
	FLOAT centerLineVelocity = flowField.getCenterLineVelocity()[i];
	// FLOAT centerLineVelocity = flowField.getVelocity().getVector(i, _parameters.geometry.sizeY / 2, _parameters.geometry.sizeZ / 2)[0];
	FLOAT downstreamRe = centerLineVelocity * x * _parameters.flow.Re;
	FLOAT delta = 0.382 * x / pow(downstreamRe, 0.2);
	return std::min(kh, 0.09 * delta);
}
