/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "parstream.H" //Gives us pout()
#include <iostream>
#include <chrono>

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "CallDoAnalysis.hpp"

// Problem specific includes:
#include "BosonStarLevel.hpp"

// Star tracking
#include "STAMR.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);


    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    STAMR st_amr;

    st_amr.m_star_tracker.initial_setup(
        sim_params.star_track_centre, sim_params.do_star_track,
        sim_params.number_of_stars, sim_params.initial_star_centres,
        sim_params.star_track_resolution, sim_params.star_track_width);
    DefaultLevelFactory<BosonStarLevel> boson_star_level_fact(st_amr,
                                                                  sim_params);
    setupAMRObject(st_amr, boson_star_level_fact);

    // Instantiate AMR interpolator for mass/GW extraction
    AMRInterpolator<Lagrange<4>> interpolator(
        st_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    st_amr.set_interpolator(&interpolator);

    // Add a scheduler to GRAMR which just calls doAnalysis on every AMRLevel
    // at time 0. It is called later in postTimeStep
    RefCountedPtr<CallDoAnalysis> call_do_analysis_ptr(new CallDoAnalysis);
    RefCountedPtr<Scheduler> scheduler_ptr(new Scheduler);
    scheduler_ptr->schedule(call_do_analysis_ptr, sim_params.max_steps);
    st_amr.schedule(scheduler_ptr);

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Engage! Run the evolution.
    st_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    st_amr.conclude();

    // Write Chombo timer report
    CH_TIMER_REPORT();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
