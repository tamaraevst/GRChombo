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
#include "RotatingBosonStarParams.hpp"
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

        pp.load("G_Newton", G_Newton, 1.0);

        pp.load("star_centre", rotating_bosonstar_params.star_centre,
                {0.5 * L, 0.5 * L, 0.5 * L});
        pp.load("initial_data_path", rotating_bosonstar_params.base_path, std::string(""));
        pp.load("BS_frequency", rotating_bosonstar_params.BS_frequency, 0.15910835770266477);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);

        //Tagging
        pp.load("star_radius", puncture_radius, 4.);
        pp.load("star_mass", puncture_mass, 4.);
        pp.load("tag_buffer", tag_buffer, 0.5);
	pp.load("tag_puncture_max_level", tag_puncture_max_level, max_level);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("mass_write_extraction",
                mass_extraction_params.write_extraction,
                false);
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
    double G_Newton;
    Real regrid_threshold_phi, regrid_threshold_chi;
    Real tag_radius_A, tag_buffer;
    Real puncture_mass, puncture_radius;
    int tag_puncture_max_level;

    RotatingBosonStar_params_t rotating_bosonstar_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;
    int flux_extraction_level; // specifies times (level) to do angmom flux extraction
};

#endif /* SIMULATIONPARAMETERS_HPP_ */