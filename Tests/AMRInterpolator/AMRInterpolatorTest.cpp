#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//General includes:
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>

#include "parstream.H" //Gives us pout()
using std::endl;
#include "AMR.H"

#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

//Problem specific includes:
#include "InterpolatorTestLevelFactory.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int
main(int argc ,char* argv[])
{
    mainSetup(argc, argv);

    //Load the parameter file and construct the SimulationParameter class
    //To add more parameters edit the SimulationParameters file.
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);
    SimulationParameters sim_params(pp);

    InterpolatorTestLevelFactory interpolator_test_level_fact(sim_params);
    AMR amr;
    setupAMRObject(amr, interpolator_test_level_fact);

    if ( status == 0 ) pout() << "GRChombo finished." << endl ;
    else pout() << "GRChombo failed with return code " << status << endl ;

    mainFinalize();
    return status ;
}
