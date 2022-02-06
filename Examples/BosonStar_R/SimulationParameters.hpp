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
#include "AngMomFluxParams.hpp"

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
        // for regridding
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("eigen", bosonstar_params.eigen, 0);
        pp.load("gridpoints",bosonstar_params.gridpoints,400000);

        pp.load("star_centre", bosonstar_params.star_centre,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);
        pp.load("BS_binary", bosonstar_params.BS_binary, false);
        pp.load("BS_BH_binary", bosonstar_params.BS_BH_binary, false);
        pp.load("BlackHoleMass", bosonstar_params.BlackHoleMass, 0.);
        pp.load("BS_rapidity", bosonstar_params.BS_rapidity, 0.0);
        pp.load("BS_separation", bosonstar_params.BS_separation, 0.0);
        pp.load("BS_impact_parameter", bosonstar_params.BS_impact_parameter, 0.0);
        pp.load("mass_ratio", bosonstar_params.mass_ratio, 1.0);
        pp.load("do_stretch", bosonstar_params.do_stretch, false);
        pp.load("n_power", bosonstar_params.n_power, 3);
        pp.load("id_choice", bosonstar_params.id_choice, 3);
        pp.load("alpha_stretch", bosonstar_params.alpha_stretch, 0.1);
        pp.load("G_Newton", bosonstar_params.Newtons_constant, 1.0);

        // Initialize values for bosonstar2_params to same as bosonstar_params
        // and then assign that ones that should differ below
        bosonstar2_params = bosonstar_params;

        // Are the two stars' profiles identical
        pp.load("identical", identical, false);

        // Boson Star 2 parameters
        if (!identical)
        {
            pp.load("central_amplitude_CSF2",
                    bosonstar2_params.central_amplitude_CSF);
            pp.load("BS_rapidity2",
                    bosonstar2_params.BS_rapidity);
        }

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Work out the minimum extraction level
        auto min_extraction_level_it = mass_extraction_params.min_extraction_level();

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);

        // Variables for outputting to plot files
        //pp.load("num_plot_vars", num_plot_vars, 0);
        //pp.load("plot_vars", plot_vars, num_plot_vars, 0);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);


        pp.load("flux_extraction_level", flux_extraction_level, 0);
        /*pp.load("flux_number_of_radii", angmomflux_params.number_radii,1);
        pp.load("flux_do", angmomflux_params.do_flux_integration,false);
        pp.load("flux_extraction_level", angmomflux_params.extraction_level,0);
        pp.load("flux_num_theta", angmomflux_params.num_theta,10);
        pp.load("flux_num_phi", angmomflux_params.num_phi,10);
        pp.load("flux_extraction_centre", angmomflux_params.centre,
                                                {0.5 * L, 0.5 * L, 0.5 * L});

        angmomflux_params.radii.resize(angmomflux_params.number_radii);
        pp.load("flux_extraction_radii", angmomflux_params.radii,
                                                angmomflux_params.number_radii);*/
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;

    // Initial data for matter and potential
    double G_Newton;
    bool identical; // whether or not the 2 boson stars have the same profile

    BosonStar_params_t bosonstar_params;
    BosonStar_params_t bosonstar2_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;


    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    // Vars for outputting in plot files
    //int num_plot_vars;
    //std::vector<int> plot_vars;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

    int flux_extraction_level; // specifies times (level) to do angmom flux extraction
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
