#ifndef _TURBULENT_FLOW_FIELD_H_
#define _TURBULENT_FLOW_FIELD_H_

#include "FlowField.h"
#include "DataStructures.h"
#include "Parameters.h"
#include "DataStructures.h"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class TurbulentFlowField : public FlowField {

    private:

        ScalarField _distance;
		    ScalarField _nitau;
        FLOAT *_centerLineBuffer;

    public:

        /** Constructor for the 2D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         */
        TurbulentFlowField (int Nx, int Ny, int Nz);

        /** Constructor for the 3D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
         */
        TurbulentFlowField ( int Nx, int Ny );

        /** Constructs a field from parameters object
         *
         * Constructs a field from a parameters object, so that it dimensionality can be defined in
         * the configuration file.
         *
         * @param parameters Parameters object with geometric information
         */
        TurbulentFlowField (const Parameters & parameters);

	    ScalarField & getTurbulentViscosity();
	    ScalarField & getWallDistance();
	    FLOAT getWallDistance(int i, int j);
	    FLOAT getWallDistance(int i, int j, int k);
        FLOAT *& getCenterLineVelocity();

        ~TurbulentFlowField ();
};

#endif
