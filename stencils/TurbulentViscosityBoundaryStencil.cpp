#include "TurbulentViscosityBoundaryStencil.h"

TurbulentViscosityBoundaryStencil::TurbulentViscosityBoundaryStencil ( const Parameters & parameters ) :
    BoundaryStencil<TurbulentFlowField> ( parameters ),
    _parameters(parameters) {}

// 2D stencils

void TurbulentViscosityBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j ){
    flowField.getTurbulentViscosity().getScalar(i, j) = flowField.getTurbulentViscosity().getScalar(i+1, j);
}


void TurbulentViscosityBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j ){
    flowField.getTurbulentViscosity().getScalar(i, j) = flowField.getTurbulentViscosity().getScalar(i-1, j);
}


void TurbulentViscosityBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){
    flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i, j+1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j+1);
}

void TurbulentViscosityBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j ){
    flowField.getTurbulentViscosity().getScalar(i, j) = -flowField.getTurbulentViscosity().getScalar(i, j-1) * _parameters.meshsize->getDy(i, j)/_parameters.meshsize->getDy(i, j-1);
}


// 3D stencils

void TurbulentViscosityBoundaryStencil::applyLeftWall ( TurbulentFlowField & flowField, int i, int j, int k ){

}


void TurbulentViscosityBoundaryStencil::applyRightWall ( TurbulentFlowField & flowField, int i, int j , int k ){

}


void TurbulentViscosityBoundaryStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){

}


void TurbulentViscosityBoundaryStencil::applyTopWall ( TurbulentFlowField & flowField, int i, int j, int k ){

}


void TurbulentViscosityBoundaryStencil::applyFrontWall ( TurbulentFlowField & flowField, int i, int j, int k ){


}


void TurbulentViscosityBoundaryStencil::applyBackWall ( TurbulentFlowField & flowField, int i, int j, int k ){

}

