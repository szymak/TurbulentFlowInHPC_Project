#include <string>
#include "VTKStencil.h"
#include "TurbulentFlowField.h"
#include "Parameters.h"

class TurbulentVTKStencil : public VTKStencil, public FieldStencil<TurbulentFlowField> {
	
	protected:
		std::stringstream _turbulentViscosityStringStream;
		std::stringstream _viscousStressXYStringStream;
		std::stringstream _viscousStressXZStringStream;
		std::stringstream _ReStressXYStringStream;
		std::stringstream _ReStressXZStringStream;
		// Store shear stress XY and XZ (mainly for channel scenario)
		FLOAT _shearStresses[2];
        FLOAT _localVelocity [ 27 * 3 ];
        FLOAT _localMeshsize [ 27 * 3 ];
	public:
		TurbulentVTKStencil(Parameters & parameters);
		void apply(TurbulentFlowField & flowField, int i, int j);
		void apply(TurbulentFlowField & flowField, int i, int j, int k);
		void write(TurbulentFlowField & flowField, int timeStep, std::string foldername);
		void writeTurbulentViscosity();
		void writeShearStresses();
		void computeShearStresses(TurbulentFlowField & flowField, int i, int j);
		void computeShearStresses(TurbulentFlowField & flowField, int i, int j, int k);
		void clearStringStreams();
		
};
