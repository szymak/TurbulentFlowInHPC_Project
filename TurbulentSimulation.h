#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#include "Simulation.h"
#include "TurbulentFlowField.h"
#include "stencils/TurbulenceFGHStencil.h"
#include "stencils/TurbulentViscosityStencil.h"
#include "stencils/TurbulentViscosityBoundaryStencil.h"
#include "stencils/TurbulentVTKStencil.h"
#include "stencils/MaxNuStencil.h"
#include "parallelManagers/PetscParallelManagerTurbulent.h"

class TurbulentSimulation : public Simulation {

	protected:
		TurbulentFlowField &_turbulentFlowField;
    	FieldIterator<TurbulentFlowField> _turbulentFghIterator;
		TurbulenceFGHStencil _turbulentFghStencil;
    	FieldIterator<TurbulentFlowField> _turbulentVtkIterator;
		TurbulentVTKStencil _turbulentVtkStencil;
		TurbulentViscosityStencil _turbulentViscosityStencil;
		FieldIterator<TurbulentFlowField> _turbulentViscosityIterator;
		MaxNuStencil _maxNuStencil;
		FieldIterator<TurbulentFlowField> _maxNuFieldIterator;
        GlobalBoundaryIterator<TurbulentFlowField> _maxNuBoundaryIterator;
        GlobalBoundaryIterator<TurbulentFlowField> _turbulentViscosityBoundaryIterator;
		TurbulentViscosityBoundaryStencil _turbulentViscosityBoundaryStencil;
		PetscParallelManagerTurbulent _parallelManagerTurbulent;

		/*
			The changes are only in the flow field and the stencils, however,
			if we change the stencils we also need to change the iterators (apparently).
			Note that _turbulentFlowField refers to the same variable as _flowField in the base class,
			but the type is TurbulentFlowField instead of just FlowField.
			So, if we need the additional functionality for turbulent flows, we access them
			through the _turbulentFlowField variable.  Otherwise, using _flowField is the same.
			All stencils operate on the same flow field.
			It does not seem to be a very elegant solution, and inheriting from Simulation does not
			seem very useful, as most of the member functions need to be reimplemented.
			Do you have any suggestions?
		*/

	public:

		TurbulentSimulation(Parameters &parameters, TurbulentFlowField &flowField);
		virtual ~TurbulentSimulation(){}
		virtual void solveTimestep();
		virtual void plotVTK(int timeStep, std::string foldername);
		virtual void initializeFlowField();



	protected:
		virtual void setTimeStep();

};

#endif // _TURBULENT_SIMULATION_H_
