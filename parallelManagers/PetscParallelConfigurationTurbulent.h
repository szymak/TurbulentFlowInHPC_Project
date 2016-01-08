#ifndef _PETSC_PARALLEL_CONFIGURATION_TURBULENT_H_
#define _PETSC_PARALLEL_CONFIGURATION_TURBULENT_H_

#include "PetscParallelConfiguration.h"

class PetscParallelConfigurationTurbulent : public PetscParallelConfiguration {

    private:
        void centerline();
        void set_centerline_flags();
    public:
        PetscParallelConfigurationTurbulent(Parameters & parameters);

};

#endif
