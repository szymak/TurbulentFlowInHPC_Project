#include "MaxNuStencil.h"
#include <algorithm>
#include <math.h>


MaxNuStencil::MaxNuStencil (const Parameters & parameters) :
    FieldStencil<TurbulentFlowField> (parameters), BoundaryStencil<TurbulentFlowField> (parameters), _Re(parameters.flow.Re) {
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
    FLOAT localMax = 2.0/((_Re + viscosity)*_factor);

    if (localMax > _maxValue){
        _maxValue = localMax;
    }
   
}

void MaxNuStencil::cellMaxValue(TurbulentFlowField & flowField, int i, int j, int k){
    FLOAT  viscosity = flowField.getTurbulentViscosity().getScalar(i, j, k);
    FLOAT localMax = 2.0/((_Re + viscosity)*_factor);

    if (localMax > _maxValue){
        _maxValue = localMax;
    }
}

void MaxNuStencil::reset () {
    _maxValue = 0;
}

void MaxNuStencil::loadFactor (FLOAT factor) {
    _factor = factor;
}

FLOAT  MaxNuStencil::getMaxValue() const{
    return _maxValue;
}
