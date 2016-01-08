#ifndef _WALL_DISTANCE_STENCIL_H_
#define _WALL_DISTANCE_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"


class WallDistanceStencil : public FieldStencil<TurbulentFlowField> {

	protected:
		
		const Parameters & _parameters;

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        WallDistanceStencil ( const Parameters & parameters );
		~WallDistanceStencil ();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        virtual void apply ( TurbulentFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        virtual void apply ( TurbulentFlowField & flowField, int i, int j, int k );

		/*
			Compute distance to the step, applicable for both 2D and 3D scenarios
			as long as the step spans the entire Z axis
		*/
		FLOAT distanceToStep(FLOAT posX, FLOAT posY);
		
};

#endif
