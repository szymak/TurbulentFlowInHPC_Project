#include "PressureBufferFillStencil.h"

PressureBufferFillStencil::PressureBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<FlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void PressureBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	this->_bufferLeftWall[j] = flowField.getPressure().getScalar(i+2, j);
}

void PressureBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j) {
	this->_bufferRightWall[j] = flowField.getPressure().getScalar(i-1, j);
}

void PressureBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	this->_bufferBottomWall[i] = flowField.getPressure().getScalar(i, j+2);
}

void PressureBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j) {
	this->_bufferTopWall[i] = flowField.getPressure().getScalar(i, j-1);
}


// 3D stencils
void PressureBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferLeftWall[flowField.getCellsZ()*j+k] = flowField.getPressure().getScalar(i+2, j, k);
}

void PressureBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferRightWall[flowField.getCellsZ()*j+k] = flowField.getPressure().getScalar(i-1, j, k);
}

void PressureBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferBottomWall[flowField.getCellsZ()*i+k] = flowField.getPressure().getScalar(i, j+2, k);
}

void PressureBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferTopWall[flowField.getCellsZ()*i+k] = flowField.getPressure().getScalar(i, j-1, k);
}

void PressureBufferFillStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferFrontWall[flowField.getCellsY()*i+j] = flowField.getPressure().getScalar(i, j, k+2);
}

void PressureBufferFillStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferBackWall[flowField.getCellsY()*i+j] = flowField.getPressure().getScalar(i, j, k-1);
}

