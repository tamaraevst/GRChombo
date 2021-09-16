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

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass, 1.0);
        //        pp.load("bh_velocity", bg_params.velocity, 0.0);
        pp.load("bh_spin", bg_params.spin, 0.0);
        pp.load("bh_center", bg_params.center, center);
        pp.load("field_amplitude", field_amplitude);

        pp.load("inner_r", inner_r, 5.0);
        pp.load("outer_r", outer_r, 500.0);

        pp.load("regrid_length", regrid_length, L);

        // Whether to do calculation of scalars' norms
        pp.load("calculate_scalar_norm", calculate_scalar_norm, false);

        // Whether to compare with analytic solution of \phi with GB term as a source (only for Schwarzschild)
        pp.load("compare_gb_analytic", compare_gb_analytic, false);

        /* Amplitudes set in front of Chern Simons and Gauss Bonnet scalars, 
        they are \gamma'(0) and \beta'(0) for the scalars respectively.
        Set them to zero if you do not want the corresponding scalar included. */
        pp.load("gamma_amplitude", gamma_amplitude, 0.0); // for Chern Simons
        pp.load("beta_amplitude", beta_amplitude, 0.0); // for Gauss Bonnet

        // directory to store data (extraction files, puncture data, constraint
        // norms)
        pp.load("data_subpath", data_path, std::string(""));
        if (!data_path.empty() && data_path.back() != '/')
            data_path += "/";
        if (output_path != "./" && !output_path.empty())
            data_path = output_path + data_path;
         
    }

    // Problem specific parameters
    double field_amplitude, sigma, regrid_length;
    int nan_check;
    double inner_r, outer_r;

    double gamma_amplitude;
    double beta_amplitude;

    bool calculate_scalar_norm;
    bool compare_gb_analytic;

    std::string data_path;


    //    std::array<double, CH_SPACEDIM> origin,
    //        dx; // location of coarsest origin and dx
    // Collection of parameters necessary for the sims
    IsotropicKerrFixedBG::params_t bg_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
