#ifndef _SIMULATION_FACTORY_H_
#define _SIMULATION_FACTORY_H_

#include "Parameters.h"
#include "Simulation.h"
#include "FlowField.h"

class SimulationFactory {

	protected:
		Simulation * _simulation;
		FlowField *_flowField;

	public:
		SimulationFactory(Parameters &parameters);
		virtual ~SimulationFactory();
		Simulation * getSimulation();

};

#endif
