#include "TurbulenceFGHStencil.h"

TurbulenceFGHStencil::TurbulenceFGHStencil ( const Parameters & parameters ) :
FieldStencil<TurbulentFlowField> ( parameters ){}


void TurbulenceFGHStencil::apply ( TurbulentFlowField & flowField,  int i, int j ){

    // Load local variables into the local array
    loadLocalVelocity2D           ( flowField  , _localVelocity          , i, j);
    loadLocalMeshsize2D           ( _parameters, _localMeshsize          , i, j);
    loadLocalTurbulentViscosity2D ( flowField  , _localTurbulentViscosity, i, j);

    FLOAT* const values = flowField.getFGH().getVector(i,j);

    values [0] = computeF2DTurbulence(_localVelocity, _localMeshsize,\
                                      _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
    values [1] = computeG2DTurbulence(_localVelocity, _localMeshsize,\
                                      _localTurbulentViscosity, _parameters, _parameters.timestep.dt);
}


void TurbulenceFGHStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ){

    const int obstacle = flowField.getFlags().getValue(i, j, k);
    FLOAT * const values = flowField.getFGH().getVector(i,j,k);

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
