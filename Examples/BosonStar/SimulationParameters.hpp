/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"

// Problem specific includes:
#include "CCZ4.hpp"
#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"

class SimulationParameters
{
  public:
    SimulationParameters(GRParmParse &pp) { readParams(pp); }

    void readParams(GRParmParse &pp)
    {
        // The automatically generated read parameters code defined in
        // SimulationParameters.inc
        auto_read_params(pp);

        // Fill in the Matter Parameters
        bosonstar_params.central_amplitude_CSF = central_amplitude_CSF;
        bosonstar_params.abs_error = abs_error;
        bosonstar_params.rel_error = rel_error;
        bosonstar_params.initial_step_size = initial_step_size;
        bosonstar_params.max_radius = max_radius;
        bosonstar_params.binary_search_tol = binary_search_tol;
        bosonstar_params.max_binary_search_iter = max_binary_search_iter;
        bosonstar_params.star_centre = star_centre;

        // Fill in the potential parameters
        potential_params.scalar_mass = scalar_mass;
        potential_params.phi4_coeff = phi4_coeff;

        // Fill in the ccz4Parameters
        ccz4_params.kappa1 = kappa1;
        ccz4_params.kappa2 = kappa2;
        ccz4_params.kappa3 = kappa3;
        ccz4_params.shift_Gamma_coeff = shift_Gamma_coeff;
        ccz4_params.shift_advec_coeff = shift_advec_coeff;
        ccz4_params.eta = eta;
        ccz4_params.lapse_power = lapse_power;
        ccz4_params.lapse_coeff = lapse_coeff;
        ccz4_params.lapse_advec_coeff = lapse_advec_coeff;
    }

// SimulationParameters.inc declares all variables and defines
// auto_read_params(GRParmParse& pp)
#include "SimulationParameters.inc"

    // Collection of parameters necessary for the CCZ4 RHS
    CCZ4::params_t ccz4_params;
    BosonStar_params_t bosonstar_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
