#include "PressureBufferReadStencil.h"

PressureBufferReadStencil::PressureBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<FlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i+1, j) = _bufferLeftWall[j];
}

void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferRightWall[j];
}

void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j+1) = _bufferBottomWall[i];
}

void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j) {
	flowField.getPressure().getScalar(i, j) = _bufferTopWall[i];
}



// 3D stencils
void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i+1, j, k) = _bufferLeftWall[flowField.getCellsZ()*j+k];
}

void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i, j, k) = _bufferRightWall[flowField.getCellsZ()*j+k];
}

void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i, j+1, k) = _bufferBottomWall[flowField.getCellsZ()*i+k];
}

void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i, j, k) = _bufferTopWall[flowField.getCellsZ()*i+k];
}

void PressureBufferReadStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i, j, k+1) = _bufferFrontWall[flowField.getCellsY()*i+j];
}

void PressureBufferReadStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {
	flowField.getPressure().getScalar(i, j, k) = _bufferBackWall[flowField.getCellsY()*i+j];
}

