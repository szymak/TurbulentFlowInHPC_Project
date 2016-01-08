#include "VelocityBufferFillStencil.h"

VelocityBufferFillStencil::VelocityBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall): BoundaryStencil<FlowField>(parameters) {
	this->_bufferLeftWall = bufferLeftWall;
	this->_bufferRightWall = bufferRightWall;
	this->_bufferTopWall = bufferTopWall;
	this->_bufferBottomWall = bufferBottomWall;
	this->_bufferFrontWall = bufferFrontWall;
	this->_bufferBackWall = bufferBackWall;
}

// 2D stencils
void VelocityBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j) {
	this->_bufferLeftWall[BUFFER_2D_POS_X(j)] = flowField.getVelocity().getVector(i+2, j, 0)[0];
	this->_bufferLeftWall[BUFFER_2D_POS_Y(j)] = flowField.getVelocity().getVector(i+2, j, 0)[1];
}

void VelocityBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j) {
	this->_bufferRightWall[BUFFER_2D_POS_X(j)] = flowField.getVelocity().getVector(i-2, j, 0)[0];
	this->_bufferRightWall[BUFFER_2D_POS_Y(j)] = flowField.getVelocity().getVector(i-1, j, 0)[1];
}

void VelocityBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j) {
	this->_bufferBottomWall[BUFFER_2D_POS_X(i)] = flowField.getVelocity().getVector(i, j+2, 0)[0];
	this->_bufferBottomWall[BUFFER_2D_POS_Y(i)] = flowField.getVelocity().getVector(i, j+2, 0)[1];
}

void VelocityBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j) {
	this->_bufferTopWall[BUFFER_2D_POS_X(i)] = flowField.getVelocity().getVector(i, j-1, 0)[0];
	this->_bufferTopWall[BUFFER_2D_POS_Y(i)] = flowField.getVelocity().getVector(i, j-2, 0)[1];
}


// 3D stencils TODO Indicies
void VelocityBufferFillStencil::applyLeftWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferLeftWall[BUFFER_3D_POS_X(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i+2, j, k)[0];
	this->_bufferLeftWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i+2, j, k)[1];
	this->_bufferLeftWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i+2, j, k)[2];
}

void VelocityBufferFillStencil::applyRightWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferRightWall[BUFFER_3D_POS_X(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i-2, j, k)[0];
	this->_bufferRightWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i-1, j, k)[1];
	this->_bufferRightWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*j+k)] = flowField.getVelocity().getVector(i-1, j, k)[2];
}

void VelocityBufferFillStencil::applyBottomWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferBottomWall[BUFFER_3D_POS_X(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j+2, k)[0];
	this->_bufferBottomWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j+2, k)[1];
	this->_bufferBottomWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j+2, k)[2];
}

void VelocityBufferFillStencil::applyTopWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferTopWall[BUFFER_3D_POS_X(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j-1, k)[0];
	this->_bufferTopWall[BUFFER_3D_POS_Y(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j-2, k)[1];
	this->_bufferTopWall[BUFFER_3D_POS_Z(flowField.getCellsZ()*i+k)] = flowField.getVelocity().getVector(i, j-1, k)[2];
}

void VelocityBufferFillStencil::applyFrontWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferFrontWall[BUFFER_3D_POS_X(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k+2)[0];
	this->_bufferFrontWall[BUFFER_3D_POS_Y(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k+2)[1];
	this->_bufferFrontWall[BUFFER_3D_POS_Z(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k+2)[2];
}

void VelocityBufferFillStencil::applyBackWall(FlowField & flowField, int i, int j, int k) {
	this->_bufferBackWall[BUFFER_3D_POS_X(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k-1)[0];
	this->_bufferBackWall[BUFFER_3D_POS_Y(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k-1)[1];
	this->_bufferBackWall[BUFFER_3D_POS_Z(flowField.getCellsY()*i+j)] = flowField.getVelocity().getVector(i, j, k-2)[2];
}

