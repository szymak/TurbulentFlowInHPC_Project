#include "SimulationFactory.h"
#include "Simulation.h"
#include "TurbulentSimulation.h"
#include "FlowField.h"
#include "TurbulentFlowField.h"

SimulationFactory::SimulationFactory(Parameters &parameters) {

	if(parameters.simulation.type == "turbulence") {
	
		_flowField = new TurbulentFlowField(parameters);
		_simulation = new TurbulentSimulation(parameters, *static_cast<TurbulentFlowField*>(_flowField));
		
	} else if(parameters.simulation.type == "dns") {
	
		_flowField = new FlowField(parameters);
		_simulation = new Simulation(parameters, *_flowField);
		
	} else {
	
		handleError(1, "Unknown simulation type! Currently supported: dns, turbulence");
      
	}
	
	if(parameters.parallel.rank == 0) {
		std::cout << "==================================================================" << std::endl;
		std::cout << "Starting simulation" << std::endl;
		std::cout << "Dimensions: " << parameters.geometry.dim << std::endl;
		std::cout << "Type: " << parameters.simulation.type << std::endl;
		std::cout << "==================================================================" << std::endl;
	}
	
}

SimulationFactory::~SimulationFactory() {
	delete _flowField;
	delete _simulation;
}

Simulation * SimulationFactory::getSimulation() {
	return _simulation;
}
