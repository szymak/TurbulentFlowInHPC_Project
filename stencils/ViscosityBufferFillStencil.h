#ifndef __VISCOSITY_BUFFER_FILL_STENCIL_H__
#define __VISCOSITY_BUFFER_FILL_STENCIL_H__

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

class ViscosityBufferFillStencil : public BoundaryStencil<TurbulentFlowField> {

	private:

	FLOAT *_bufferLeftWall;
	FLOAT *_bufferRightWall;
	FLOAT *_bufferTopWall;
	FLOAT *_bufferBottomWall;
	FLOAT *_bufferFrontWall;
	FLOAT *_bufferBackWall;

	public:

	ViscosityBufferFillStencil (const Parameters & parameters, FLOAT *bufferLeftWall, FLOAT *bufferRightWall, FLOAT *bufferTopWall, FLOAT *bufferBottomWall, FLOAT *bufferFrontWall, FLOAT *bufferBackWall);

	void applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j);

	void applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j);

	void applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j);

	void applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j);

	void applyLeftWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

	void applyRightWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

	void applyBottomWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

	void applyTopWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

	void applyFrontWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

	void applyBackWall(TurbulentFlowField & turbulentfFowField, int i, int j, int k);

};

#endif
