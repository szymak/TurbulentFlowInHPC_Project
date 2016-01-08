#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "Configuration.h"
#include "SimulationFactory.h"
#include "Simulation.h"
#include "TurbulentSimulation.h"
#include "parallelManagers/PetscParallelConfiguration.h"
#include "parallelManagers/PetscParallelConfigurationTurbulent.h"
#include "MeshsizeFactory.h"
#include <time.h>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "TurbulentFlowField.h"  // TODO move this line to the new TurbulentSimulation

int main (int argc, char *argv[]) {

    // Parallelization related. Initialize and identify
    // ---------------------------------------------------
    int rank;   // This processor's identifier
    int nproc;  // Number of processors in the group
    clock_t start, end; // Variables to measure elapsed time
    PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::cout << "Rank: " << rank << ", Nproc: " << nproc << std::endl;
    //----------------------------------------------------


    // read configuration and store information in parameters object
    Configuration configuration(argv[1]);
    Parameters parameters;
    configuration.loadParameters(parameters);
    PetscParallelConfiguration * parallelConfiguration;
    if (parameters.simulation.type=="dns")
      parallelConfiguration = new  PetscParallelConfiguration (parameters);
    else if (parameters.simulation.type=="turbulence")
      parallelConfiguration = new  PetscParallelConfigurationTurbulent (parameters);


    MeshsizeFactory::getInstance().initMeshsize(parameters);

    SimulationFactory simulationFactory(parameters);
    Simulation *simulation = simulationFactory.getSimulation();

    #ifdef DEBUG
    std::cout << "Processor " << parameters.parallel.rank << " with index ";
    std::cout << parameters.parallel.indices[0] << ",";
    std::cout << parameters.parallel.indices[1] << ",";
    std::cout << parameters.parallel.indices[2];
    std::cout <<    " is computing the size of its subdomain and obtains ";
    std::cout << parameters.parallel.localSize[0] << ", ";
    std::cout << parameters.parallel.localSize[1] << " and ";
    std::cout << parameters.parallel.localSize[2];
    std::cout << ". Left neighbour: " << parameters.parallel.leftNb;
    std::cout << ", right neighbour: " << parameters.parallel.rightNb;
    std::cout << std::endl;
    std::cout << "Min. meshsizes: " << parameters.meshsize->getDxMin() << ", " << parameters.meshsize->getDyMin() << ", " << parameters.meshsize->getDzMin() << std::endl;
    #endif

    // call initialization of simulation (initialize flow field)
    if(simulation == NULL){ handleError(1, "simulation==NULL!"); }
    simulation->initializeFlowField();
    //flowField->getFlags().show();

    /* Create a folder for the VTK output */ // TODO: Shouldn't this part be executed only by the rank0 process?
	time_t rawtime = time(NULL);
	struct tm * timeinfo = localtime(&rawtime);
	char buffer[80];
	strftime (buffer, sizeof(buffer),"%d.%m.%y-%H:%M:%S", timeinfo);
	std::stringstream foldername;
	foldername << parameters.vtk.prefix << "_" << buffer;
	struct stat info;
	if( stat( foldername.str().c_str(), &info ) != 0 ) {
		mkdir(foldername.str().c_str(), 0777);
		/* Copy the currently used config file into the VTK output directory */
		const char* fileName = "/config.xml";
		std::ifstream source(argv[1], std::ios::binary);
		std::ofstream dest(foldername.str().append(fileName).c_str(), std::ios::binary);
    		dest << source.rdbuf();
    		source.close();
    		dest.close();
	}

    FLOAT time = 0.0;
    FLOAT timeStdOut=parameters.stdOut.interval;
    FLOAT timeVTKOut=parameters.vtk.interval;
    int timeSteps = 0;

    // TODO WS1: plot initial state
    simulation->plotVTK(timeSteps, foldername.str());
    timeVTKOut += parameters.vtk.interval;

    // START: measure program time
    if (rank == 0) {   
    	start = clock(); 
	}
    // time loop
    while (time < parameters.simulation.finalTime){

      simulation->solveTimestep();

      time += parameters.timestep.dt;

      // std-out: terminal info
      if ( (rank==0) && (timeStdOut <= time) ){
          std::cout << "Current time: " << time << "\ttimestep: " <<
                        parameters.timestep.dt << std::endl;
          timeStdOut += parameters.stdOut.interval;
      }

      timeSteps++;

      // TODO WS1: trigger VTK output
      if(timeVTKOut <= time) {
	      simulation->plotVTK(timeSteps, foldername.str());
	      timeVTKOut += parameters.vtk.interval;
      }

    }

    // STOP: finish measurement
    if (rank == 0) {
    	end = clock();
    	std::cout << "TOTAL TIME: " << (float)(end-start) / CLOCKS_PER_SEC << std::endl;
    }
    // TODO WS1: plot final output
    simulation->plotVTK(timeSteps, foldername.str());

    if (rank == 0) {
        /* Copy the output file in the results folder */
        const char* fileName = "/MYJOB.out";
        std::ifstream source("MYJOB.out", std::ios::binary);
        std::ofstream dest(foldername.str().append(fileName).c_str(), std::ios::binary);
        dest << source.rdbuf();
        source.close();
        dest.close();
    }

    PetscFinalize();

    delete parallelConfiguration;
}
