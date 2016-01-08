#include "MaxNuStencil.h"
#include <algorithm>
#include <math.h>


MaxNuStencil::MaxNuStencil (const Parameters & parameters) :
    FieldStencil<TurbulentFlowField> (parameters), BoundaryStencil<TurbulentFlowField> (parameters) {
    reset();
}

void MaxNuStencil::apply (TurbulentFlowField & flowField, int i, int j){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::apply (TurbulentFlowField & flowField, int i, int j, int k){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxNuStencil::applyLeftWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyRightWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyBottomWall ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyTopWall    ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyFrontWall  ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxNuStencil::applyBackWall   ( TurbulentFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}


void MaxNuStencil::cellMaxValue(TurbulentFlowField & flowField, int i, int j){
    FLOAT  viscosity = flowField.getTurbulentViscosity().getScalar(i, j);
   
    if (viscosity > _maxValue){
        _maxValue = viscosity;
    }
   
}

void MaxNuStencil::cellMaxValue(TurbulentFlowField & flowField, int i, int j, int k){
    FLOAT  viscosity= flowField.getTurbulentViscosity().getScalar(i, j, k);
   
    if (viscosity > _maxValue){
        _maxValue = viscosity;
    }
}

void MaxNuStencil::reset () {
    _maxValue = 0;
}

FLOAT  MaxNuStencil::getMaxValue() const{
    return _maxValue;
}
