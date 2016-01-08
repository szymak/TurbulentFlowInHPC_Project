#include "ViscosityBufferFillStencil.h"

ViscosityBufferFillStencil::ViscosityBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<TurbulentFlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j) {
	this->_bufferLeftWall[j] = turbulentfFowField.getTurbulentViscosity().getScalar(i+2, j);
}

void ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j) {
	this->_bufferRightWall[j] = turbulentfFowField.getTurbulentViscosity().getScalar(i-1, j);
}

void ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j) {
	this->_bufferBottomWall[i] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j+2);
}

void ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j) {
	this->_bufferTopWall[i] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j-1);
}


// 3D stencils
void ViscosityBufferFillStencil::applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferLeftWall[turbulentfFowField.getCellsZ()*j+k] = turbulentfFowField.getTurbulentViscosity().getScalar(i+2, j, k);
}

void ViscosityBufferFillStencil::applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferRightWall[turbulentfFowField.getCellsZ()*j+k] = turbulentfFowField.getTurbulentViscosity().getScalar(i-1, j, k);
}

void ViscosityBufferFillStencil::applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferBottomWall[turbulentfFowField.getCellsZ()*i+k] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j+2, k);
}

void ViscosityBufferFillStencil::applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferTopWall[turbulentfFowField.getCellsZ()*i+k] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j-1, k);
}

void ViscosityBufferFillStencil::applyFrontWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferFrontWall[turbulentfFowField.getCellsY()*i+j] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j, k+2);
}

void ViscosityBufferFillStencil::applyBackWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k) {
	this->_bufferBackWall[turbulentfFowField.getCellsY()*i+j] = turbulentfFowField.getTurbulentViscosity().getScalar(i, j, k-1);
}
