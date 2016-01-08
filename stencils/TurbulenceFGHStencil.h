#ifndef _TURBULENCE_STENCIL_FGH_H_
#define _TURBULENCE_STENCIL_FGH_H_

#include "../TurbulentFlowField.h"
#include "../Stencil.h"
#include "../Parameters.h"
#include "StencilFunctions.h"
#include "Definitions.h"


class TurbulenceFGHStencil : public FieldStencil<TurbulentFlowField>
{

    private:
      //Local variable that will be used to approximate derivatives.
        FLOAT _localVelocity           [ 27 * 3 ];
        FLOAT _localMeshsize           [ 27 * 3 ];
        FLOAT _localTurbulentViscosity [ 27 * 3 ];      //    TODO shoud I change it to 27 instead of 27*3

    public:
        //constructor
        TurbulenceFGHStencil ( const Parameters & parameters );
        //apply methods
        void apply ( TurbulentFlowField & flowField, int i, int j );
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );
};


#endif
