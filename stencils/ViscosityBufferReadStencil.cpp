#include "ViscosityBufferReadStencil.h"

ViscosityBufferReadStencil::ViscosityBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField & turbulentFlowField, int i, int j) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i+1, j) = _bufferLeftWall[j];
}

void ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField & turbulentFlowField, int i, int j) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j) = _bufferRightWall[j];
}

void ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField & turbulentFlowField, int i, int j) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j+1) = _bufferBottomWall[i];
}

void ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField & turbulentFlowField, int i, int j) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j) = _bufferTopWall[i];
}



// 3D stencils
void ViscosityBufferReadStencil::applyLeftWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i+1, j, k) = _bufferLeftWall[turbulentFlowField.getCellsZ()*j+k];
}

void ViscosityBufferReadStencil::applyRightWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j, k) = _bufferRightWall[turbulentFlowField.getCellsZ()*j+k];
}

void ViscosityBufferReadStencil::applyBottomWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j+1, k) = _bufferBottomWall[turbulentFlowField.getCellsZ()*i+k];
}

void ViscosityBufferReadStencil::applyTopWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j, k) = _bufferTopWall[turbulentFlowField.getCellsZ()*i+k];
}

void ViscosityBufferReadStencil::applyFrontWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j, k+1) = _bufferFrontWall[turbulentFlowField.getCellsY()*i+j];
}

void ViscosityBufferReadStencil::applyBackWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k) {
	turbulentFlowField.getTurbulentViscosity().getScalar(i, j, k) = _bufferBackWall[turbulentFlowField.getCellsY()*i+j];
}
