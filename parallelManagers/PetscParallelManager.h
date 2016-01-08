#ifndef __PETSC_PARALLEL_MANAGER_H__
#define __PETSC_PARALLEL_MANAGER_H__

#include "../FlowField.h"
#include "../Parameters.h"
#include "../Iterators.h"
#include "../stencils/PressureBufferFillStencil.h"
#include "../stencils/PressureBufferReadStencil.h"
#include "../stencils/VelocityBufferFillStencil.h"
#include "../stencils/VelocityBufferReadStencil.h"

class PetscParallelManager {

	protected:
	
	int _cellsX;
	int _cellsY;
	int _cellsZ; // Set to 1 for 2D case
	
	// Number of elements in the parallel boundaries
	int _cellsLeftRight;
	int _cellsTopBottom;
	int _cellsFrontBack;
	
	Parameters &_parameters;
	FlowField &_flowField;
	
	// Buffers
	FLOAT *_pressureSendBufferLeftWall;
	FLOAT *_pressureRecvBufferLeftWall;
	FLOAT *_pressureSendBufferRightWall;
	FLOAT *_pressureRecvBufferRightWall;
	FLOAT *_pressureSendBufferTopWall;
	FLOAT *_pressureRecvBufferTopWall;
	FLOAT *_pressureSendBufferBottomWall;
	FLOAT *_pressureRecvBufferBottomWall;
	FLOAT *_pressureSendBufferFrontWall;
	FLOAT *_pressureRecvBufferFrontWall;
	FLOAT *_pressureSendBufferBackWall;
	FLOAT *_pressureRecvBufferBackWall;
	FLOAT *_velocitySendBufferLeftWall;
	FLOAT *_velocityRecvBufferLeftWall;
	FLOAT *_velocitySendBufferRightWall;
	FLOAT *_velocityRecvBufferRightWall;
	FLOAT *_velocitySendBufferTopWall;
	FLOAT *_velocityRecvBufferTopWall;
	FLOAT *_velocitySendBufferBottomWall;
	FLOAT *_velocityRecvBufferBottomWall;
	FLOAT *_velocitySendBufferFrontWall;
	FLOAT *_velocityRecvBufferFrontWall;
	FLOAT *_velocitySendBufferBackWall;
	FLOAT *_velocityRecvBufferBackWall;
	
	// Stencils
	PressureBufferFillStencil _pressureBufferFillStencil;
	PressureBufferReadStencil _pressureBufferReadStencil;
	VelocityBufferFillStencil _velocityBufferFillStencil;
	VelocityBufferReadStencil _velocityBufferReadStencil;
	
	// Iterators
	ParallelBoundaryIterator<FlowField> _pressureBufferFillIterator;
	ParallelBoundaryIterator<FlowField> _pressureBufferReadIterator;
	ParallelBoundaryIterator<FlowField> _velocityBufferFillIterator;
	ParallelBoundaryIterator<FlowField> _velocityBufferReadIterator;

	public:
	
	PetscParallelManager(FlowField &flowField, Parameters &parameters);
	~PetscParallelManager();
	void communicatePressure();
	void communicateVelocity();
	void sendReceive(FLOAT *sendBuffer, int sendTo, FLOAT *receiveBuffer, int receiveFrom, int size);
	

};

#endif
