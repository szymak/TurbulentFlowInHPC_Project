#include "CenterLineVelocityBufferFillStencil.h"

CenterLineVelocityBufferFillStencil::CenterLineVelocityBufferFillStencil (const Parameters & parameters, FLOAT *centerLineBuffer): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_centerLineBuffer = centerLineBuffer;
}

// 2D stencils
void CenterLineVelocityBufferFillStencil::applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j) {}

void CenterLineVelocityBufferFillStencil::applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j) {}

void CenterLineVelocityBufferFillStencil::applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j) {

	this->_centerLineBuffer[i] =  turbulentfFowField.getVelocity().getVector(i, j+_parameters.parallel.local_center_line_index[1])[0];

}

void CenterLineVelocityBufferFillStencil::applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j) {}


// 3D stencils
void CenterLineVelocityBufferFillStencil::applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {}

void CenterLineVelocityBufferFillStencil::applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {}

void CenterLineVelocityBufferFillStencil::applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {

	this->_centerLineBuffer[i] =  turbulentfFowField.getVelocity().getVector(i, j+_parameters.parallel.local_center_line_index[1],k+_parameters.parallel.local_center_line_index[2])[0];

}

void CenterLineVelocityBufferFillStencil::applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {}

void CenterLineVelocityBufferFillStencil::applyFrontWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {}

void CenterLineVelocityBufferFillStencil::applyBackWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {}
