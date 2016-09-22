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

//Problem specific includes:
#include "CCZ4LevelFactory.hpp"

//TODO (MK): There is a lot of clutter still in this file ... get rid of everything that's not necessary
//and comment on the rest (so that we don't carry around useless code forever

#include <fenv.h>
#if defined(__i386__) && defined(__SSE__)
#include <xmmintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

int
runGRChombo (ParmParse& pp)
{
    //The line below selects the problem that is simulated
    //(every problem must have a child of AMRLevel which is produced
    //in a child of AMRLevelFactory)
    CCZ4LevelFactory ccz4_level_fact(pp);
    AMR amr;
    setupAMRObject(amr, ccz4_level_fact);

    Real stop_time;
    pp.get("stop_time", stop_time);
    int max_steps;
    pp.get("max_steps", max_steps);

    amr.run(stop_time, max_steps);

    amr.conclude();

    return 0 ;
}

int
main(int argc ,char* argv[])
{
    mainSetup(argc, argv);

    char* in_file = argv[1];
 //TODO: this whole ParmParse business is annoying. Read the parameters then pass SimulationParamters around
    ParmParse  pp(argc-2,argv+2,NULL,in_file);
    int status = runGRChombo(pp);

    if ( status == 0 ) pout() << "GRChombo finished." << endl ;
    else pout() << "GRChombo failed with return code " << status << endl ;

    mainFinalize();
    return status ;
}
