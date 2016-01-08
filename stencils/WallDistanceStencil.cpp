#include <algorithm>
#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "VTKStencil.h"
#include "Iterators.h"
#include "WallDistanceStencil.h"

WallDistanceStencil::WallDistanceStencil ( const Parameters & parameters ) : 
	FieldStencil<TurbulentFlowField> ( parameters ),
	_parameters(parameters)
{
}

WallDistanceStencil::~WallDistanceStencil () {
}

void WallDistanceStencil::apply ( TurbulentFlowField & flowField, int i, int j ) {
	
	FLOAT dx = _parameters.meshsize->getDx(i, j);
	FLOAT dy = _parameters.meshsize->getDy(i, j);
	FLOAT posX = fabs(_parameters.meshsize->getPosX(i, j) + 0.5*dx);
	FLOAT posY = fabs(_parameters.meshsize->getPosY(i, j) + 0.5*dy);
	
	// Min distance to top or bottom wall
	FLOAT minDistance = std::min(posY, _parameters.geometry.lengthY-posY);
	
	if(this->_parameters.simulation.scenario == "channel" && _parameters.bfStep.xRatio > 0 && _parameters.bfStep.yRatio > 0) {
		
		// Consider the backward facing step
		minDistance = std::min(minDistance, distanceToStep(posX, posY));
		
	} else if(this->_parameters.simulation.scenario == "cavity") {
		
		// Consider also the lateral walls
		FLOAT minX = std::min(posX, _parameters.geometry.lengthX-posX);
		minDistance = std::min(minDistance, minX);
	}

	flowField.getWallDistance().getScalar(i, j) = minDistance;
}

void WallDistanceStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {

	FLOAT dx = _parameters.meshsize->getDx(i, j, k);
	FLOAT dy = _parameters.meshsize->getDy(i, j, k);
	FLOAT dz = _parameters.meshsize->getDz(i, j, k);
	FLOAT posX = fabs(_parameters.meshsize->getPosX(i, j, k) + 0.5*dx);
	FLOAT posY = fabs(_parameters.meshsize->getPosY(i, j, k) + 0.5*dy);
	FLOAT posZ = fabs(_parameters.meshsize->getPosZ(i, j, k) + 0.5*dz);
	FLOAT minY = std::min(posY, _parameters.geometry.lengthY-posY);
	FLOAT minZ = std::min(posZ, _parameters.geometry.lengthZ-posZ);
	
	// Min distance to top, bottom, front or back wall
	FLOAT minDistance = std::min(minY, minZ);
	
	if(this->_parameters.simulation.scenario == "channel" && _parameters.bfStep.xRatio > 0 && _parameters.bfStep.yRatio > 0) {
		
		// Consider the backward facing step
		minDistance = std::min(minDistance, distanceToStep(posX, posY));
		
	} else if(this->_parameters.simulation.scenario == "cavity") {

		// Consider also the lateral walls
		FLOAT minX = std::min(posX, _parameters.geometry.lengthX-posX);
		minDistance = std::min(minDistance, minX);
	}

	flowField.getWallDistance().getScalar(i, j, k) = minDistance;
}

FLOAT WallDistanceStencil::distanceToStep(FLOAT posX, FLOAT posY) {
	FLOAT stepWidth = _parameters.bfStep.xRatio * _parameters.geometry.lengthX;
	FLOAT stepHeight = _parameters.bfStep.yRatio * _parameters.geometry.lengthY;
	FLOAT distanceToStep;
	if(posX < stepWidth) {
		// above the step
		distanceToStep = posY - stepHeight;
	} else if(posX > stepWidth && posY < stepHeight) {
		// right of the step
		distanceToStep = posX - stepWidth;
	} else {
		// diagonal to the step (compute distance to corner)
		distanceToStep = sqrt((posX - stepWidth)*(posX - stepWidth) + (posY - stepHeight)*(posY - stepHeight));
	}
	return distanceToStep;
}

