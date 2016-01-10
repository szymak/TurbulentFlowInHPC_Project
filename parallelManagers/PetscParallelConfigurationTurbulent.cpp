#include "PetscParallelConfigurationTurbulent.h"

PetscParallelConfigurationTurbulent::PetscParallelConfigurationTurbulent(Parameters & parameters):
	PetscParallelConfiguration(parameters)
{
	MPI_Comm_split(PETSC_COMM_WORLD, _parameters.parallel.firstCorner[0], _parameters.parallel.rank, &_parameters.parallel.planeComm);
	MPI_Comm_rank(_parameters.parallel.planeComm, &_parameters.parallel.plane_rank);

	set_centerline_flags();
	centerline();
}

void PetscParallelConfigurationTurbulent::centerline(){

		int numPlanesX = _parameters.parallel.numProcessors[0];
		int myPlaneId = _parameters.parallel.indices[0];
		int * localRootOfPlane = new int[numPlanesX]();
		int * globalRootOfPlane = new int[numPlanesX]();

		localRootOfPlane[myPlaneId] = _parameters.parallel.centerlineFlag * _parameters.parallel.plane_rank;

		MPI_Allreduce(localRootOfPlane, globalRootOfPlane, numPlanesX, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		_parameters.parallel.plane_root = globalRootOfPlane[myPlaneId];

		int globalYCenterIdx = (_parameters.geometry.sizeY + 3) / 2;
		int localYCenterIdx = globalYCenterIdx - _parameters.parallel.firstCorner[1];

		_parameters.parallel.local_center_line_index[1] = localYCenterIdx;

		int globalZCenterIdx = (_parameters.geometry.sizeZ + 3) / 2;

		int localZCenterIdx = globalZCenterIdx - _parameters.parallel.firstCorner[2];
		_parameters.parallel.local_center_line_index[2] = localZCenterIdx;

		delete [] localRootOfPlane;
		delete [] globalRootOfPlane;

}

void PetscParallelConfigurationTurbulent::set_centerline_flags()
{
  int center_processor [3]={};

  _parameters.parallel.centerlineFlag=0;

  _parameters.parallel.local_center_line_index[0]=0;
  _parameters.parallel.local_center_line_index[1]=0;
  _parameters.parallel.local_center_line_index[2]=0;

  //determine the center processor in X direction
  center_processor[0]=_parameters.parallel.indices[0];

  //determine the center processor in Y direction
  if (_parameters.parallel.indices[0] < _parameters.bfStep.xRatio * _parameters.parallel.numProcessors[0] ) {
    center_processor[1] = floor ( ( _parameters.parallel.numProcessors[1] * ( 1 + _parameters.bfStep.yRatio ) ) / 2 );
  }
  else
  {
    center_processor[1] = floor ( ( _parameters.geometry.lengthY / 2 ) / ( _parameters.geometry.lengthY / _parameters.parallel.numProcessors[1] ) );
  }
  _parameters.parallel.centerProcessor[1]=center_processor[1];

  //determine the center processor in Z direction
  center_processor[2] = floor ( _parameters.parallel.numProcessors[2] / 2.0 );
  _parameters.parallel.centerProcessor[2]=center_processor[2];

  //determie the local index of the centerline in the center processor
  if ( ( center_processor[1]==_parameters.parallel.indices[1] ) && ( center_processor[2]==_parameters.parallel.indices[2] ) )
  {
    _parameters.parallel.centerlineFlag=1;

    //in Y direction
    if ( ( (_parameters.parallel.numProcessors[1] + lrint( _parameters.bfStep.yRatio * _parameters.parallel.numProcessors[1] ) ) % 2 ) == 1 )
    {
      _parameters.parallel.local_center_line_index[1]=_parameters.parallel.localSize[1]/2;
    }
    else
    {
      _parameters.parallel.local_center_line_index[1]=0;
    }

    //in Z direction
    if ( ( _parameters.parallel.numProcessors[2] % 2 ) == 1 )
    {
      _parameters.parallel.local_center_line_index[2]=_parameters.parallel.localSize[2]/2;
    }
    else
    {
      _parameters.parallel.local_center_line_index[2]=0;
    }

  }

}
