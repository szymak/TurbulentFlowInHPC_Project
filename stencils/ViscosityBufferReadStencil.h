#ifndef __VISCOSITY_BUFFER_READ_STENCIL_H__
#define __VISCOSITY_BUFFER_READ_STENCIL_H__

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

class ViscosityBufferReadStencil : public BoundaryStencil<TurbulentFlowField> {

	private:

	FLOAT *_bufferLeftWall;
	FLOAT *_bufferRightWall;
	FLOAT *_bufferTopWall;
	FLOAT *_bufferBottomWall;
	FLOAT *_bufferFrontWall;
	FLOAT *_bufferBackWall;

	public:

	ViscosityBufferReadStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall);

	void applyLeftWall(TurbulentFlowField & turbulentFlowField, int i, int j);

	void applyRightWall(TurbulentFlowField & turbulentFlowField, int i, int j);

	void applyBottomWall(TurbulentFlowField & turbulentFlowField, int i, int j);

	void applyTopWall(TurbulentFlowField & turbulentFlowField, int i, int j);

	void applyLeftWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

	void applyRightWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

	void applyBottomWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

	void applyTopWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

	void applyFrontWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

	void applyBackWall(TurbulentFlowField & turbulentFlowField, int i, int j, int k);

};


#endif
