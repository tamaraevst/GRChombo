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
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include "InterpolationQuery.hpp"
#include "UserVariables.hpp"

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

    //Setup the AMRInterpolator
    int num_points = 2;

    double* chi_ptr = new double[num_points];
    double* phi_ptr = new double[num_points];
    double* interp_x = new double[num_points];
    double* interp_y = new double[num_points];
    double* interp_z = new double[num_points];

    double L;
    pp.get("L", L);
    double extract_center[3] = {L/2,L/2,L/2};
    double extract_radius = L/4;

    for (int iPhi=0; iPhi<num_points; ++iPhi)
    {
        double extract_angle = iPhi*M_PI/num_points;
        interp_x[iPhi] = extract_center[0] + extract_radius*cos(extract_angle);
        interp_y[iPhi] = extract_center[1] + extract_radius*sin(extract_angle);
        interp_z[iPhi] = extract_center[2];
    }

    InterpolationQuery query(num_points);
    query
        .setCoords(0, interp_x)
        .setCoords(1, interp_y)
        .setCoords(2, interp_z)
        .addComp(c_chi, chi_ptr)
        .addComp(c_phi, phi_ptr);

    double dx_scalar = GRAMRLevel::gr_cast(amr.getAMRLevels()[0])->get_dx(); //coarsest grid spacing
    Array<double, 3> dx;
    Array<double, CH_SPACEDIM> origin;
    for (int idir=0; idir<CH_SPACEDIM; ++idir)
    {
        dx[idir] = dx_scalar;
        origin[idir] = dx_scalar/2;
    }

    AMRInterpolator<Lagrange<4> > interpolator(amr, origin, dx, 2);
    interpolator.interp(query);

    int status = 0;

    if ( status == 0 ) pout() << "GRChombo finished." << endl ;
    else pout() << "GRChombo failed with return code " << status << endl ;

    mainFinalize();
    return status;
}
