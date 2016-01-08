#ifndef __CENTER_LINE_VELOCITY_BUFFER_FILL_STENCIL_H__
#define __CENTER_LINE_VELOCITY_BUFFER_FILL_STENCIL_H__

#include "../Stencil.h"
#include "../TurbulentFlowField.h"

class CenterLineVelocityBufferFillStencil : public BoundaryStencil<TurbulentFlowField> {

	private:

	FLOAT *_centerLineBuffer;

	public:

	CenterLineVelocityBufferFillStencil (const Parameters & parameters, FLOAT *_centerLineBuffer);

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
