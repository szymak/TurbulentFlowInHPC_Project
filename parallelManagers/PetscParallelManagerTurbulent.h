#ifndef __PETSC_PARALLEL_MANAGER_TURBULENT_H__
#define __PETSC_PARALLEL_MANAGER_TURBULENT_H__

#include "PetscParallelManager.h"
#include "../stencils/ViscosityBufferFillStencil.h"
#include "../stencils/ViscosityBufferReadStencil.h"
#include "../stencils/CenterLineVelocityBufferFillStencil.h"
//#include "../TurbulentFlowField.h"

class PetscParallelManagerTurbulent: public PetscParallelManager {

  private:
  TurbulentFlowField &_turbulentFlowField;

  FLOAT *_viscositySendBufferLeftWall;
  FLOAT *_viscosityRecvBufferLeftWall;
  FLOAT *_viscositySendBufferRightWall;
  FLOAT *_viscosityRecvBufferRightWall;
  FLOAT *_viscositySendBufferTopWall;
  FLOAT *_viscosityRecvBufferTopWall;
  FLOAT *_viscositySendBufferBottomWall;
  FLOAT *_viscosityRecvBufferBottomWall;
  FLOAT *_viscositySendBufferFrontWall;
  FLOAT *_viscosityRecvBufferFrontWall;
  FLOAT *_viscositySendBufferBackWall;
  FLOAT *_viscosityRecvBufferBackWall;

  FLOAT *_centerLineBuffer;

  ViscosityBufferFillStencil _viscosityBufferFillStencil;
	ViscosityBufferReadStencil _viscosityBufferReadStencil;
  CenterLineVelocityBufferFillStencil _centerLineVelocityFillStencil;

  ParallelBoundaryIterator<TurbulentFlowField> _viscosityBufferFillIterator;
  ParallelBoundaryIterator<TurbulentFlowField> _viscosityBufferReadIterator;
  ParallelBoundaryIterator<TurbulentFlowField> _centerLineVelocityFillIterator;


  public:

  PetscParallelManagerTurbulent(TurbulentFlowField &flowField, Parameters &parameters);
  ~PetscParallelManagerTurbulent();
  void communicateViscosity();
  void communicateCenterLineVelocity();
};

#endif
