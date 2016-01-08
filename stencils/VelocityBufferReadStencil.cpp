#include "VelocityBufferReadStencil.h"

VelocityBufferReadStencil::VelocityBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<FlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void VelocityBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferLeftWall[BUFFER_2D_POS_X(j)];
	flowField.getVelocity().getVector(i+1, j, 0)[1] = _bufferLeftWall[BUFFER_2D_POS_Y(j)];
}

void VelocityBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferRightWall[BUFFER_2D_POS_X(j)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferRightWall[BUFFER_2D_POS_Y(j)];
}

void VelocityBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j+1, 0)[0] = _bufferBottomWall[BUFFER_2D_POS_X(i)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferBottomWall[BUFFER_2D_POS_Y(i)];
}

void VelocityBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j) {
	flowField.getVelocity().getVector(i, j, 0)[0] = _bufferTopWall[BUFFER_2D_POS_X(i)];
	flowField.getVelocity().getVector(i, j, 0)[1] = _bufferTopWall[BUFFER_2D_POS_Y(i)];
}


// 3D stencils
void VelocityBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j, k)[0] = _bufferLeftWall[BUFFER_3D_POS_X(flowField.getCellsZ()*j+k)];
	flowField.getVelocity().getVector(i+1, j, k)[1] = _bufferLeftWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*j+k)];
	flowField.getVelocity().getVector(i+1, j, k)[2] = _bufferLeftWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*j+k)];
}

void VelocityBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j, k)[0] = _bufferRightWall[BUFFER_3D_POS_X(flowField.getCellsZ()*j+k)];
	flowField.getVelocity().getVector(i, j, k)[1] = _bufferRightWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*j+k)];
	flowField.getVelocity().getVector(i, j, k)[2] = _bufferRightWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*j+k)];
}

void VelocityBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j+1, k)[0] = _bufferBottomWall[BUFFER_3D_POS_X(flowField.getCellsZ()*i+k)];
	flowField.getVelocity().getVector(i, j, k)[1] = _bufferBottomWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*i+k)];
	flowField.getVelocity().getVector(i, j+1, k)[2] = _bufferBottomWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*i+k)];
}

void VelocityBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j, k)[0] = _bufferTopWall[BUFFER_3D_POS_X(flowField.getCellsZ()*i+k)];
	flowField.getVelocity().getVector(i, j, k)[1] = _bufferTopWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*i+k)];
	flowField.getVelocity().getVector(i, j, k)[2] = _bufferTopWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*i+k)];
}

void VelocityBufferReadStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j, k+1)[0] = _bufferFrontWall[BUFFER_3D_POS_X(flowField.getCellsY()*i+j)];
	flowField.getVelocity().getVector(i, j, k+1)[1] = _bufferFrontWall[BUFFER_3D_POS_Y(flowField.getCellsY()*i+j)];
	flowField.getVelocity().getVector(i, j, k)[2] = _bufferFrontWall[BUFFER_3D_POS_Z(flowField.getCellsY()*i+j)];
}

void VelocityBufferReadStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {
	flowField.getVelocity().getVector(i, j, k)[0] = _bufferBackWall[BUFFER_3D_POS_X(flowField.getCellsY()*i+j)];
	flowField.getVelocity().getVector(i, j, k)[1] = _bufferBackWall[BUFFER_3D_POS_Y(flowField.getCellsY()*i+j)];
	flowField.getVelocity().getVector(i, j, k)[2] = _bufferBackWall[BUFFER_3D_POS_Z(flowField.getCellsY()*i+j)];
}

