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
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "ScalarField.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction
        pp.load("scalar_amplitude", initial_params.amplitude, 0.1);
        pp.load("scalar_width", initial_params.width, 1.0);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass, 1.0);
        pp.load("kerr_spin", kerr_params.spin, 0.0);
        pp.load("kerr_center", kerr_params.center, center);

        // whether to do extraction
        pp.load("activate_extraction", activate_extraction, false);

        // Whether to do calculation of scalars' norms
        pp.load("calculate_scalar_norm", calculate_scalar_norm, false);

        // Whether to compare with analytic solution of \phi with GB term as a source (only for Schwarzschild)
        pp.load("compare_gb_analytic", compare_gb_analytic, false);

        //Whether to evolve just matter counterpart
        pp.load("evolve_just_matter", evolve_just_matter, false);

        /* Amplitudes set in front of Chern Simons and Gauss Bonnet scalars, 
        they are \gamma'(0) and \beta'(0) for the scalars respectively.
        Set them to zero if you do not want the corresponding scalar included. */
        pp.load("gamma_amplitude", gamma_amplitude, 0.0); // for Chern Simons
        pp.load("beta_amplitude", beta_amplitude, 0.0); // for Gauss Bonnet
    }

    void check_params()
    {
        warn_parameter("scalar_width", initial_params.width,
                       initial_params.width < 0.5 * L,
                       "is greater than half the domain size");
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    InitialScalarData::params_t initial_params;
    KerrBH::params_t kerr_params;

    //Parameters for modified scalar field equation 
    double gamma_amplitude;
    double beta_amplitude;
    bool calculate_scalar_norm;
    bool compare_gb_analytic;
    bool evolve_just_matter;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
