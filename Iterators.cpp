template<class FlowField>
FieldIterator<FlowField>::FieldIterator (FlowField & flowField, const Parameters& parameters, FieldStencil<FlowField> & stencil,
                               int lowOffset, int highOffset):
    Iterator<FlowField>(flowField,parameters), _stencil(stencil), _lowOffset(lowOffset), _highOffset(highOffset){}


template<class FlowField>
void FieldIterator<FlowField>::iterate (){

    const int cellsX = Iterator<FlowField>::_flowField.getCellsX();
    const int cellsY = Iterator<FlowField>::_flowField.getCellsY();
    const int cellsZ = Iterator<FlowField>::_flowField.getCellsZ();
    // The index k can be used for the 2D and 3D cases.

    if (Iterator<FlowField>::_parameters.geometry.dim == 2){

        // Loop without lower boundaries. These will be dealt with by the global boundary stencils
        // or by the subdomain boundary iterators.
        for (int j = 1 + _lowOffset; j < cellsY - 1 + _highOffset; j++){
            for (int i = 1 + _lowOffset; i < cellsX - 1 + _highOffset; i++){
                _stencil.apply ( Iterator<FlowField>::_flowField, i, j );
            }
        }
    }

    if (Iterator<FlowField>::_parameters.geometry.dim == 3){
        for (int k = 1 + _lowOffset; k < cellsZ - 1 + _highOffset; k++){
            for (int j = 1 + _lowOffset; j < cellsY - 1 + _highOffset; j++){
                for (int i = 1 + _lowOffset; i < cellsX - 1 + _highOffset; i++){
                    _stencil.apply ( Iterator<FlowField>::_flowField, i, j, k );
                }
            }
        }
    }
}


template<class FlowField>
OMPIterator<FlowField>::OMPIterator (FlowField & flowField, const Parameters& parameters, FieldStencil<FlowField> & stencil,
                               int lowOffset, int highOffset):
    Iterator<FlowField>(flowField,parameters), _stencil(stencil), _lowOffset(lowOffset), _highOffset(highOffset){}


template<class FlowField>
void OMPIterator<FlowField>::iterate (){

    const int cellsX = Iterator<FlowField>::_flowField.getCellsX();
    const int cellsY = Iterator<FlowField>::_flowField.getCellsY();
    const int cellsZ = Iterator<FlowField>::_flowField.getCellsZ();
    // The index k can be used for the 2D and 3D cases.

    if (Iterator<FlowField>::_parameters.geometry.dim == 2){

        // Loop without lower boundaries. These will be dealt with by the global boundary stencils
        // or by the subdomain boundary iterators.
        for (int j = 1 + _lowOffset; j < cellsY - 1 + _highOffset; j++){
            for (int i = 1 + _lowOffset; i < cellsX - 1 + _highOffset; i++){
                _stencil.apply ( Iterator<FlowField>::_flowField, i, j );
            }
        }
    }

    if (Iterator<FlowField>::_parameters.geometry.dim == 3){	
	        for (int k = 1 + _lowOffset; k < cellsZ - 1 + _highOffset; k++){
		    	//black    
		        #pragma omp parallel for schedule (static)
		            for (int j = 1 + _lowOffset; j < cellsY - 1 + _highOffset; j++){
		                for (int i = 1 + _lowOffset + (j + (k%2))%2; i < cellsX - 1 + _highOffset; i+=2){
                            apply ( Iterator<FlowField>::_flowField, i, j, k );
		                }
		            }
			  	//red
				#pragma omp parallel for schedule (static)
		            for (int j = 1 + _lowOffset; j < cellsY - 1 + _highOffset; j++){
		                for (int i = 2 + _lowOffset - (j + (k%2))%2; i < cellsX - 1 + _highOffset; i+=2){
                            apply ( Iterator<FlowField>::_flowField, i, j, k );
		                }
		            }
	        }	
    }
    
}


template<class FlowField>
void OMPIterator<FlowField>::apply (TurbulentFlowField & flowField, int i, int j, int k ){

    const int obstacle = flowField.getFlags().getValue(i, j, k);
    FLOAT * const values = flowField.getFGH().getVector(i,j,k);

    FLOAT _localVelocity           [ 27 * 3];
    FLOAT _localMeshsize           [ 27 * 3];
    FLOAT _localTurbulentViscosity [ 27 * 3];      //    TODO shoud I change it to 27 instead of 27*3


    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid

        loadLocalVelocity3D(  flowField, _localVelocity,                      i, j, k);
        loadLocalMeshsize3D(_parameters, _localMeshsize,                      i, j, k);
        loadLocalTurbulentViscosity3D ( flowField , _localTurbulentViscosity, i, j, k);

        if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
            values [0] = computeF3DTurbulence(_localVelocity, _localMeshsize,\
                                              _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            values [1] = computeG3DTurbulence(_localVelocity, _localMeshsize,\
                                              _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            values [2] = computeH3DTurbulence(_localVelocity, _localMeshsize,\
                                              _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
        }
    }

}

template<class FlowField>
GlobalBoundaryIterator<FlowField>::GlobalBoundaryIterator(FlowField & flowField,
                                               const Parameters & parameters,
                                               BoundaryStencil<FlowField> & stencil,
                                               int lowOffset, int highOffset):
    Iterator<FlowField> (flowField,parameters),
    _lowOffset(lowOffset), _highOffset(highOffset),
    _leftWallStencil(stencil), _rightWallStencil(stencil),
    _bottomWallStencil(stencil), _topWallStencil(stencil),
    _frontWallStencil(stencil), _backWallStencil(stencil){}


template<class FlowField>
GlobalBoundaryIterator<FlowField>::GlobalBoundaryIterator(FlowField & flowField,
                                               const Parameters & parameters,
                                               BoundaryStencil<FlowField> & leftWallStencil,
                                               BoundaryStencil<FlowField> & rightWallStencil,
                                               BoundaryStencil<FlowField> & bottomWallStencil,
                                               BoundaryStencil<FlowField> & topWallStencil,
                                               int lowOffset, int highOffset):
    Iterator<FlowField> (flowField,parameters),
    _lowOffset(lowOffset), _highOffset(highOffset),
    _leftWallStencil(leftWallStencil), _rightWallStencil(rightWallStencil),
    _bottomWallStencil(bottomWallStencil), _topWallStencil(topWallStencil),
    // This is plain bad, but it will work. The references had to be initialized somehow
    _frontWallStencil(leftWallStencil), _backWallStencil(leftWallStencil)
{

    if (Iterator<FlowField>::_parameters.geometry.dim == 3){
        handleError(1, "Trying to use 2D constructor for a 3D field");
    }

}


template<class FlowField>
GlobalBoundaryIterator<FlowField>::GlobalBoundaryIterator(FlowField & flowField,
                                               const Parameters & parameters,
                                               BoundaryStencil<FlowField> & leftWallStencil,
                                               BoundaryStencil<FlowField> & rightWallStencil,
                                               BoundaryStencil<FlowField> & bottomWallStencil,
                                               BoundaryStencil<FlowField> & topWallStencil,
                                               BoundaryStencil<FlowField> & frontWallStencil,
                                               BoundaryStencil<FlowField> & backWallStencil,
                                               int lowOffset, int highOffset):
    Iterator<FlowField> (flowField,parameters),
    _lowOffset(lowOffset), _highOffset(highOffset),
    _leftWallStencil(leftWallStencil), _rightWallStencil(rightWallStencil),
    _bottomWallStencil(bottomWallStencil), _topWallStencil(topWallStencil),
    _frontWallStencil(frontWallStencil), _backWallStencil(backWallStencil){}



template<class FlowField>
void GlobalBoundaryIterator<FlowField>::iterate () {

    if (Iterator<FlowField>::_parameters.geometry.dim == 2){

        if (Iterator<FlowField>::_parameters.parallel.leftNb < 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY() + _highOffset; j++) {
                _leftWallStencil.applyLeftWall (Iterator<FlowField>::_flowField, _lowOffset, j);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.rightNb < 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY() + _highOffset; j++) {
                _rightWallStencil.applyRightWall (Iterator<FlowField>::_flowField, Iterator<FlowField>::_flowField.getCellsX()+_highOffset-1,j);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.bottomNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX() + _highOffset; i++) {
                _bottomWallStencil.applyBottomWall (Iterator<FlowField>::_flowField, i, _lowOffset);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.topNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX() + _highOffset; i++) {
                _topWallStencil.applyTopWall (Iterator<FlowField>::_flowField, i, Iterator<FlowField>::_flowField.getCellsY()+_highOffset-1);
            }
        }
    }

    if (Iterator<FlowField>::_parameters.geometry.dim == 3){

        if (Iterator<FlowField>::_parameters.parallel.leftNb < 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _leftWallStencil.applyLeftWall   ( Iterator<FlowField>::_flowField, _lowOffset, j, k );
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.rightNb < 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _rightWallStencil.applyRightWall (Iterator<FlowField>::_flowField, Iterator<FlowField>::_flowField.getCellsX()+_highOffset-1, j, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.bottomNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _bottomWallStencil.applyBottomWall (Iterator<FlowField>::_flowField, i, _lowOffset, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.topNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _topWallStencil.applyTopWall (Iterator<FlowField>::_flowField, i, Iterator<FlowField>::_flowField.getCellsY()+_highOffset-1, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.frontNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                    _frontWallStencil.applyFrontWall (Iterator<FlowField>::_flowField, i, j, _lowOffset);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.backNb < 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                    _backWallStencil.applyBackWall (Iterator<FlowField>::_flowField, i, j, Iterator<FlowField>::_flowField.getCellsZ()+_highOffset-1);
                }
            }
        }
    }
}


template <class FlowField>
ParallelBoundaryIterator<FlowField>::ParallelBoundaryIterator (FlowField & flowField,
                                                               const Parameters & parameters,
                                                               BoundaryStencil<FlowField> & stencil,
                                                               int lowOffset, int highOffset):
    Iterator<FlowField>(flowField,parameters), _stencil(stencil),
    _lowOffset(lowOffset), _highOffset(highOffset){}


template<class FlowField>
void ParallelBoundaryIterator<FlowField>::iterate () {

    if (Iterator<FlowField>::_parameters.geometry.dim == 2){

        if (Iterator<FlowField>::_parameters.parallel.leftNb >= 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY() + _highOffset; j++) {
                _stencil.applyLeftWall (Iterator<FlowField>::_flowField, _lowOffset, j);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.rightNb >= 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY() + _highOffset; j++) {
                _stencil.applyRightWall (Iterator<FlowField>::_flowField, Iterator<FlowField>::_flowField.getCellsX()+_highOffset-1,j);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.bottomNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX() + _highOffset; i++) {
                _stencil.applyBottomWall (Iterator<FlowField>::_flowField, i, _lowOffset);
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.topNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX() + _highOffset; i++) {
                _stencil.applyTopWall (Iterator<FlowField>::_flowField, i, Iterator<FlowField>::_flowField.getCellsY()+_highOffset-1);
            }
        }
    }

    if (Iterator<FlowField>::_parameters.geometry.dim == 3){

        if (Iterator<FlowField>::_parameters.parallel.leftNb >= 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _stencil.applyLeftWall   ( Iterator<FlowField>::_flowField, _lowOffset, j, k );
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.rightNb >= 0){
            for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _stencil.applyRightWall (Iterator<FlowField>::_flowField, Iterator<FlowField>::_flowField.getCellsX()+_highOffset-1, j, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.bottomNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _stencil.applyBottomWall (Iterator<FlowField>::_flowField, i, _lowOffset, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.topNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int k = _lowOffset; k < Iterator<FlowField>::_flowField.getCellsZ()+_highOffset; k++) {
                    _stencil.applyTopWall (Iterator<FlowField>::_flowField, i, Iterator<FlowField>::_flowField.getCellsY()+_highOffset-1, k);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.frontNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                    _stencil.applyFrontWall (Iterator<FlowField>::_flowField, i, j, _lowOffset);
                }
            }
        }

        if (Iterator<FlowField>::_parameters.parallel.backNb >= 0){
            for (int i = _lowOffset; i < Iterator<FlowField>::_flowField.getCellsX()+_highOffset; i++) {
                for (int j = _lowOffset; j < Iterator<FlowField>::_flowField.getCellsY()+_highOffset; j++) {
                    _stencil.applyBackWall (Iterator<FlowField>::_flowField, i, j, Iterator<FlowField>::_flowField.getCellsZ()+_highOffset-1);
                }
            }
        }
    }
}
