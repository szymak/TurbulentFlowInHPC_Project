#ifndef _TURBVISC_STENCIL_H_
#define _TURBVISC_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbulentFlowField.h"
#include "MixingLengthModel.h"



class TurbulentViscosityStencil : public FieldStencil<TurbulentFlowField> {

	
    private:

        // A local velocity variable that will be used to approximate derivatives. Size matches 3D
        // case, but can be used for 2D as well.
        FLOAT _localVelocity [ 27 * 3 ];
        // local meshsize
        FLOAT _localMeshsize [ 27 * 3 ];
        MixingLengthModel * _mixingLength;
        
    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        TurbulentViscosityStencil ( const Parameters & parameters );
		~TurbulentViscosityStencil ();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( TurbulentFlowField & flowField, int i, int j, int k );

   
};

#endif
