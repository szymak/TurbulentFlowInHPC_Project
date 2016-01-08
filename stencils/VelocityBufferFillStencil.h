#ifndef __VELOCITY_BUFFER_FILL_STENCIL_H__
#define __VELOCITY_BUFFER_FILL_STENCIL_H__

#include "../Stencil.h"
#include "../FlowField.h"
#define BUFFER_2D_POS_X(i) 2*(i)
#define BUFFER_2D_POS_Y(i) 2*(i)+1
#define BUFFER_3D_POS_X(i) 3*(i)
#define BUFFER_3D_POS_Y(i) 3*(i)+1
#define BUFFER_3D_POS_Z(i) 3*(i)+2

class VelocityBufferFillStencil : public BoundaryStencil<FlowField> {

	private:

	FLOAT *_bufferLeftWall;
	FLOAT *_bufferRightWall;
	FLOAT *_bufferTopWall;
	FLOAT *_bufferBottomWall;
	FLOAT *_bufferFrontWall;
	FLOAT *_bufferBackWall;

	public:

	VelocityBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall);

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
