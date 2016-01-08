#include "PetscParallelManagerTurbulent.h"
#include "PetscParallelConfiguration.h"

PetscParallelManagerTurbulent::PetscParallelManagerTurbulent(TurbulentFlowField &flowField, Parameters &parameters) :
  //parents' constructor
  PetscParallelManager(flowField, parameters),
  _turbulentFlowField(flowField),

  //Buffers
  _viscositySendBufferLeftWall		(new FLOAT[_cellsLeftRight]),
  _viscosityRecvBufferLeftWall		(new FLOAT[_cellsLeftRight]),
  _viscositySendBufferRightWall	  (new FLOAT[_cellsLeftRight]),
  _viscosityRecvBufferRightWall	  (new FLOAT[_cellsLeftRight]),
  _viscositySendBufferTopWall		  (new FLOAT[_cellsTopBottom]),
  _viscosityRecvBufferTopWall	  	(new FLOAT[_cellsTopBottom]),
  _viscositySendBufferBottomWall	(new FLOAT[_cellsTopBottom]),
  _viscosityRecvBufferBottomWall	(new FLOAT[_cellsTopBottom]),
  _viscositySendBufferFrontWall	  (parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscosityRecvBufferFrontWall	  (parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscositySendBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _viscosityRecvBufferBackWall		(parameters.geometry.dim == 2 ? NULL : new FLOAT[_cellsFrontBack]),
  _centerLineBuffer             	(new FLOAT[_cellsX]),

  //Stencils
  _viscosityBufferFillStencil(parameters, _viscositySendBufferLeftWall, _viscositySendBufferRightWall, _viscositySendBufferTopWall, _viscositySendBufferBottomWall, _viscositySendBufferFrontWall, _viscositySendBufferBackWall),
  _viscosityBufferReadStencil(parameters, _viscosityRecvBufferLeftWall, _viscosityRecvBufferRightWall, _viscosityRecvBufferTopWall, _viscosityRecvBufferBottomWall, _viscosityRecvBufferFrontWall, _viscosityRecvBufferBackWall),

  _centerLineVelocityFillStencil(parameters, _centerLineBuffer),


  //Iterators
  _viscosityBufferFillIterator(flowField, parameters, _viscosityBufferFillStencil),
  _viscosityBufferReadIterator(flowField, parameters, _viscosityBufferReadStencil),

  _centerLineVelocityFillIterator(flowField, parameters, _centerLineVelocityFillStencil,0,0)
{

}

PetscParallelManagerTurbulent::~PetscParallelManagerTurbulent(){
  delete [] _viscositySendBufferLeftWall;
  delete [] _viscosityRecvBufferLeftWall;
  delete [] _viscositySendBufferRightWall;
  delete [] _viscosityRecvBufferRightWall;
  delete [] _viscositySendBufferTopWall;
  delete [] _viscosityRecvBufferTopWall;
  delete [] _viscositySendBufferBottomWall;
  delete [] _viscositySendBufferFrontWall;
  delete [] _viscosityRecvBufferFrontWall;
  delete [] _viscositySendBufferBackWall;
  delete [] _viscosityRecvBufferBackWall;
  delete [] _centerLineBuffer;
}

void PetscParallelManagerTurbulent::communicateViscosity() {

	_viscosityBufferFillIterator.iterate();

	// Left to right & Right to left
	sendReceive(_viscositySendBufferRightWall, _parameters.parallel.rightNb, _viscosityRecvBufferLeftWall, _parameters.parallel.leftNb, _cellsLeftRight);
	sendReceive(_viscositySendBufferLeftWall, _parameters.parallel.leftNb, _viscosityRecvBufferRightWall, _parameters.parallel.rightNb, _cellsLeftRight);
	// Top to bottom & Bottom to top
	sendReceive(_viscositySendBufferTopWall, _parameters.parallel.topNb, _viscosityRecvBufferBottomWall, _parameters.parallel.bottomNb, _cellsTopBottom);
	sendReceive(_viscositySendBufferBottomWall, _parameters.parallel.bottomNb, _viscosityRecvBufferTopWall, _parameters.parallel.topNb, _cellsTopBottom);
	// Front to back & Back to front
	sendReceive(_viscositySendBufferFrontWall, _parameters.parallel.frontNb, _viscosityRecvBufferBackWall, _parameters.parallel.backNb, _cellsFrontBack);
	sendReceive(_viscositySendBufferBackWall, _parameters.parallel.backNb, _viscosityRecvBufferFrontWall, _parameters.parallel.frontNb, _cellsFrontBack);

	_viscosityBufferReadIterator.iterate();

}


void PetscParallelManagerTurbulent::communicateCenterLineVelocity() {
  // buffer fill . iterate  for centerline
/*  if (_parameters.parallel.centerlineFlag)
  {
      _centerLineVelocityFillIterator.iterate();
  }

  */
  if (_parameters.parallel.centerlineFlag && _parameters.geometry.dim == 2)
    {
  	for(int i = 0; i < _flowField.getCellsX(); i++) {
  		_centerLineBuffer[i] = _flowField.getVelocity().getVector(i, _parameters.parallel.local_center_line_index[1])[0];
  	}
    } else if (_parameters.parallel.centerlineFlag && _parameters.geometry.dim == 3) {
  	for(int i = 0; i < _flowField.getCellsX(); i++) {
  		_centerLineBuffer[i] = _flowField.getVelocity().getVector(i, _parameters.parallel.local_center_line_index[1],  _parameters.parallel.local_center_line_index[2])[0];
  	}
    }

  // communicate the center line velocity
  MPI_Bcast(
      _centerLineBuffer,
      _parameters.parallel.localSize[0]+2,
      (sizeof(FLOAT) == sizeof(float) ? MPI_FLOAT : MPI_DOUBLE),
      _parameters.parallel.plane_root,
      _parameters.parallel.planeComm);

  _turbulentFlowField.getCenterLineVelocity()=_centerLineBuffer;

}
