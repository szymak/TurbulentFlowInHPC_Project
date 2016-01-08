#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "StencilFunctions.h"
#include "../TurbulentFlowField.h"
#include "TurbulentViscosityStencil.h"
#include "Iterators.h"
#include <iomanip>
#include <mpi.h>



TurbulentViscosityStencil::TurbulentViscosityStencil ( const Parameters & parameters ) :
	FieldStencil<TurbulentFlowField> ( parameters )
{
	// Choose mixing length model from the beginning,
	// to avoid checking inside every call to "apply"
	if(_parameters.simulation.scenario == "channel" && (_parameters.bfStep.xRatio <= 0 || _parameters.bfStep.yRatio <= 0)) {
		if (_parameters.turbulence.boundary_layer_equation == "laminar") {
			_mixingLength = new LmLaminarFlatPlate(_parameters);
		} else if(_parameters.turbulence.boundary_layer_equation == "turbulent") {
			_mixingLength = new LmTurbulentFlatPlate(_parameters);
		} else if(_parameters.turbulence.boundary_layer_equation == "kh") {
			_mixingLength = new LmKh(_parameters);
		}
	} else {
		_mixingLength = new LmKh(_parameters);
		if(_parameters.turbulence.boundary_layer_equation != "kh") {
			std::cout << "!!! Warning: Using Lm = Kh as mixing length" << std::endl;
		}
	}
}

TurbulentViscosityStencil::~TurbulentViscosityStencil () {
	delete _mixingLength;	
}

void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j ) {

	FLOAT lm = _mixingLength->at(flowField, i, j);
	loadLocalVelocity2D(  flowField, _localVelocity, i, j);
	loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);
	FLOAT SdotS = computeSdotS2D(_localVelocity, _localMeshsize);

	flowField.getTurbulentViscosity().getScalar(i, j) = lm * lm * sqrt(2 * SdotS);

}

void TurbulentViscosityStencil::apply ( TurbulentFlowField & flowField, int i, int j, int k ) {

	FLOAT lm = _mixingLength->at(flowField, i, j, k);
	loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
	loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);
	FLOAT SdotS = computeSdotS3D(_localVelocity, _localMeshsize);

	flowField.getTurbulentViscosity().getScalar(i, j, k) = lm * lm * sqrt(2 * SdotS);

}
