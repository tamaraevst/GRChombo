/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "IsotropicKerrFixedBG.hpp"

class SimulationParameters : public ChomboParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : ChomboParameters(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);
        pp.load("regrid_length", regrid_length, L);

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass, 1.0);
        //        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_spin", bg_params.spin, 0.0);
        pp.load("bh_center", bg_params.center, center);
        pp.load("field_amplitude", field_amplitude);
        pp.load("inner_r", inner_r, 5.0);
        pp.load("outer_r", outer_r, 500.0);
    }

    // Problem specific parameters
    double field_amplitude, regrid_length;
    double sigma, excision_width;
    int nan_check;
    double inner_r, outer_r;
    //    std::array<double, CH_SPACEDIM> origin,
    //        dx; // location of coarsest origin and dx
    // Collection of parameters necessary for the sims
    IsotropicKerrFixedBG::params_t bg_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
