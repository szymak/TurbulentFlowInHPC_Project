#ifndef _DTW_STENCIL_H_
#define _DTW_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>



class DistanceToWallStencil : public FieldStencil<TurbulentFlowField> {


    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        VTKStencil ( const Parameters & parameters );
		~VTKStencil ();

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( FlowField & flowField, int timeStep, std::string foldername );

        std::string getFilename( int timeStep, std::string foldername );
        void writeFileHeader();
		void writeCellDataHeader ( FlowField & flowField );
        void writeGrid ( FlowField & flowField );
        void writePressure ( FlowField & flowField );
        void writeVelocity ( FlowField & flowField );

};

#endif
