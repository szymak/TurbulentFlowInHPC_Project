#ifndef _VTK_STENCIL_H_
#define _VTK_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>

/** TODO WS1: Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class VTKStencil : public FieldStencil<FlowField> {

	protected:
	
		std::ofstream *_outputFile;
		std::stringstream _velocityStringStream;
		std::stringstream _pressureStringStream;

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
        virtual void apply ( FlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        virtual void apply ( FlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        virtual void write ( FlowField & flowField, int timeStep, std::string foldername );

        virtual std::string getFilename( int timeStep, std::string foldername );
        virtual void writeFileHeader();
		virtual void writeCellDataHeader ( FlowField & flowField );
        virtual void writeGrid ( FlowField & flowField );
        virtual void writePressure ( );
        virtual void writeVelocity ( );
        virtual void clearStringStreams();

};

#endif
