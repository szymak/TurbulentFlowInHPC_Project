#ifndef __PRESSURE_BUFFER_READ_STENCIL_H__
#define __PRESSURE_BUFFER_READ_STENCIL_H__

#include "../Stencil.h"
#include "../FlowField.h"

class PressureBufferReadStencil : public BoundaryStencil<FlowField> {

	private:
	
	FLOAT *_bufferLeftWall;
	FLOAT *_bufferRightWall;
	FLOAT *_bufferTopWall;
	FLOAT *_bufferBottomWall;
	FLOAT *_bufferFrontWall;
	FLOAT *_bufferBackWall;
	
	public:

	PressureBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall);
	
	void applyLeftWall(FlowField & flowField, int i, int j);

	void applyRightWall(FlowField & flowField, int i, int j);

	void applyBottomWall(FlowField & flowField, int i, int j);

	void applyTopWall(FlowField & flowField, int i, int j);

	void applyLeftWall(FlowField & flowField, int i, int j, int k);

	void applyRightWall(FlowField & flowField, int i, int j, int k);

	void applyBottomWall(FlowField & flowField, int i, int j, int k);

	void applyTopWall(FlowField & flowField, int i, int j, int k);

	void applyFrontWall(FlowField & flowField, int i, int j, int k);

	void applyBackWall(FlowField & flowField, int i, int j, int k);
	
};


#endif
