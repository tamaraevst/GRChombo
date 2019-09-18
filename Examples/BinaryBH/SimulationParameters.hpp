/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "BoostedBH.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp)
    {
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("centerA", bh1_params.center);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("centerB", bh2_params.center);
        pp.load("momentumB", bh2_params.momentum);
        pp.load("activate_extraction", activate_extraction, 0);
        pp.load("track_punctures", track_punctures, 0);

        puncture_coords.resize(2 * CH_SPACEDIM);
        FOR1(idir)
        {
            puncture_coords[0 + idir] = bh1_params.center[idir];
            puncture_coords[CH_SPACEDIM + idir] = bh2_params.center[idir];
        }
    }

    // Initial data
    int activate_extraction, track_punctures;
    std::vector<double> puncture_coords;
    Real massA, massB;
    std::array<double, CH_SPACEDIM> centerA, centerB;
    std::array<double, CH_SPACEDIM> momentumA, momentumB;

    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
