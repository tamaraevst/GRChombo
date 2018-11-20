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
#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // Read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // Read in all the parameters from the parameter file
        auto_read_params(pp);

        // Fill in the boson star parameters struct
        bosonstar_params.central_amplitude_CSF = central_amplitude_CSF;
        bosonstar_params.abs_error = abs_error;
        bosonstar_params.rel_error = rel_error;
        bosonstar_params.initial_step_size = initial_step_size;
        bosonstar_params.max_radius = max_radius;
        bosonstar_params.binary_search_tol = binary_search_tol;
        bosonstar_params.max_binary_search_iter = max_binary_search_iter;
        bosonstar_params.star_centre = star_centre;

        // Fill in the potential parameters struct
        potential_params.scalar_mass = scalar_mass;
        potential_params.phi4_coeff = phi4_coeff;
    }

    void auto_read_params(GRParmParse& pp)
    {
        pp.load("regrid_threshold_K", regrid_threshold_K);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);

        // Newton's constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Boson Star initial data params
        pp.load("central_amplitude_CSF", central_amplitude_CSF); //this will be divided by sqrt(4pi)
        pp.load("abs_error", abs_error, 1.0e-14);
        pp.load("rel_error", rel_error, 1.0e-14);
        pp.load("initial_step_size", initial_step_size, 0.015625);
        pp.load("max_radius", max_radius, 256.);
        pp.load("binary_search_tol", binary_search_tol, 1.0e-15);
        pp.load("max_binary_search_iter", max_binary_search_iter, 1e3);
        pp.load("star_centre", star_centre, {0.5*L, 0.5*L, 0.5*L});

        // Potential params
        pp.load("scalar_mass", scalar_mass, 1.0);
        pp.load("phi4_coeff", phi4_coeff, 0.0);
    }

    // Tagging thresholds
    Real regrid_threshold_K, regrid_threshold_phi;

    // Initial data for matter and potential
    double G_Newton;
    double central_amplitude_CSF, abs_error, rel_error, initial_step_size;
    double max_radius, binary_search_tol, max_binary_search_iter;
    double scalar_mass, phi4_coeff;
    std::array<double, CH_SPACEDIM> star_centre;

    // Parameter structs
    BosonStar_params_t bosonstar_params;
    Potential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
