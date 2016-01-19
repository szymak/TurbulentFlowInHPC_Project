#include "TurbulentSimulation.h"
#include "stencils/WallDistanceStencil.h"

TurbulentSimulation::TurbulentSimulation(Parameters &parameters, TurbulentFlowField &flowField):
	Simulation(parameters, flowField),
	_turbulentFlowField(flowField),
	_turbulentFghStencil(parameters),
	_turbulentFghIterator(_turbulentFlowField, parameters, _turbulentFghStencil),
	_turbulentViscosityStencil(parameters),
	_turbulentViscosityIterator(_turbulentFlowField, parameters, _turbulentViscosityStencil, 1, 0),
	_turbulentVtkStencil(parameters),
	_turbulentVtkIterator(_turbulentFlowField, parameters, _turbulentVtkStencil, 1, 0),
	_maxNuStencil(parameters),
    _maxNuFieldIterator(_turbulentFlowField,parameters,_maxNuStencil),
    _maxNuBoundaryIterator(_turbulentFlowField,parameters,_maxNuStencil),
	_turbulentViscosityBoundaryStencil(parameters),
	_turbulentViscosityBoundaryIterator(_turbulentFlowField, parameters, _turbulentViscosityBoundaryStencil, 1, 0),
	_parallelManagerTurbulent(_turbulentFlowField, parameters)
{
}

// TODO: Changes for turbulent simulation not implemented yet!
void TurbulentSimulation::solveTimestep(){
	SimpleTimer timer;
    FLOAT totalTime; // Variable to measure elapsed time for each running process

	// determine and set max. timestep which is allowed in this simulation
	setTimeStep();

	// compute fgh
	if (_parameters.parallel.rank == 0) {timer.start();}
	_turbulentFghIterator.iterate();
    if (_parameters.parallel.rank == 0) {
        totalTime = timer.getTimeAndRestart();
        //std::cout << "FGH time step: " << totalTime << std::endl;
        iterator_times[FGH] += totalTime;
    }

	// set global boundary values
	_wallFGHIterator.iterate();

	// compute the right hand side
    if (_parameters.parallel.rank == 0) {totalTime = timer.getTimeAndRestart();}
	_rhsIterator.iterate();
    if (_parameters.parallel.rank == 0) {
        totalTime = timer.getTimeAndRestart();
        //std::cout << "RHS time step: " << totalTime << std::endl;
        iterator_times[RHS] += totalTime;
    }

	// solve for pressure
	_solver.solve();
	// TODO WS2: communicate pressure values
	_parallelManagerTurbulent.communicatePressure();

    if (_parameters.parallel.rank == 0) {totalTime = timer.getTimeAndRestart();}
	// compute velocity
	_velocityIterator.iterate();
    if (_parameters.parallel.rank == 0) {
        totalTime = timer.getTimeAndRestart();
        //std::cout << "VELO time step: " << totalTime << std::endl;
        iterator_times[VELO] += totalTime;
    }

    if (_parameters.parallel.rank == 0) {totalTime = timer.getTimeAndRestart();}
	// set obstacle boundaries
	_obstacleIterator.iterate();
    if (_parameters.parallel.rank == 0) {
        totalTime = timer.getTimeAndRestart();
        //std::cout << "OBST time step: " << totalTime << std::endl;
        iterator_times[OBST] += totalTime;
    }

	// TODO WS2: communicate velocity values
	_parallelManagerTurbulent.communicateVelocity();
	// Iterate for velocities on the boundary
	_wallVelocityIterator.iterate();

	_parallelManagerTurbulent.communicateCenterLineVelocity();

    if (_parameters.parallel.rank == 0) {totalTime = timer.getTimeAndRestart();}
	_turbulentViscosityIterator.iterate();
    if (_parameters.parallel.rank == 0) {
        totalTime = timer.getTimeAndRestart();
        //std::cout << "VISC time step: " << totalTime << std::endl;
        iterator_times[VISC] += totalTime;
    }

	_parallelManagerTurbulent.communicateViscosity();
	_turbulentViscosityBoundaryIterator.iterate();
}

// TODO: Changes for turbulent simulation not implemented yet!
void TurbulentSimulation::setTimeStep(){

	const FLOAT cinematicviscosity=1.0/_parameters.flow.Re;
	FLOAT localMin, globalMin;
	assertion(_parameters.geometry.dim == 2 || _parameters.geometry.dim == 3);
	FLOAT factor = 1.0/(_parameters.meshsize->getDxMin() * _parameters.meshsize->getDxMin()) +
		         1.0/(_parameters.meshsize->getDyMin() * _parameters.meshsize->getDyMin());

	// determine maximum velocity
	_maxUStencil.reset();
	_maxUFieldIterator.iterate();
	_maxUBoundaryIterator.iterate();

	_maxNuStencil.reset();
	_maxNuFieldIterator.iterate();
	_maxNuBoundaryIterator.iterate();

	if (_parameters.geometry.dim == 3) {
	factor += 1.0/(_parameters.meshsize->getDzMin() * _parameters.meshsize->getDzMin());
	_parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[2];
	} else {
	_parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[0];
	}

	localMin = std::min(_parameters.timestep.dt,
		                            std::min(std::min(1.0/(cinematicviscosity+_maxNuStencil.getMaxValue())/(2*factor),
		                            1.0 / _maxUStencil.getMaxValues()[0]),
		                            1.0 / _maxUStencil.getMaxValues()[1]));

	// Here, we select the type of operation before compiling. This allows to use the correct
	// data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
	// machines.

	globalMin = MY_FLOAT_MAX;
	MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

	_parameters.timestep.dt = globalMin;
	_parameters.timestep.dt *= _parameters.timestep.tau;

}

void TurbulentSimulation::plotVTK(int timeStep, std::string foldername) {
	_turbulentVtkIterator.iterate();
	_turbulentVtkStencil.write(this->_turbulentFlowField, timeStep, foldername);
}

void TurbulentSimulation::initializeFlowField() {
	Simulation::initializeFlowField();
	WallDistanceStencil wds(_parameters);
	FieldIterator<TurbulentFlowField> it(_turbulentFlowField, _parameters, wds);
	it.iterate();
};
