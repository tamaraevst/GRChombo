/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "BosonStarLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
#include "GammaCalculator.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"
#include "ConstraintViolations.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "BosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

// For mass extraction
#include "ADMMass.hpp"
#include "Density.hpp"
#include "EMTensor.hpp"
#include "MomFluxCalc.hpp"
#include "SourceIntPreconditioner.hpp"
//#include "MassExtraction.hpp"

// For GW extraction
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// For Noether Charge calculation
#include "SmallDataIO.hpp"
#include "NoetherCharge.hpp"

// For Star Tracking
#include "GaussianFitTracking.hpp"

// For Ang Mom Integrating
#include "AngMomFlux.hpp"

// Things to do at each advance step, after the RK4 is calculated
void BosonStarLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void BosonStarLevel::initialData()
{
    CH_TIME("BosonStarLevel::initialData");
    if (m_verbosity)
        pout() << "BosonStarLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    BosonStar boson_star(m_p.bosonstar_params, m_p.potential_params,
                         m_p.G_Newton, m_dx, m_verbosity);


    // the max radius the code might need to calculate out to is L*sqrt(3)
    boson_star.compute_1d_solution(2.*m_p.L);

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(GammaCalculator(m_dx),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                   disable_simd());

    fillAllGhosts();
}

// Things to do before outputting a checkpoint file
void BosonStarLevel::preCheckpointLevel()
{
    CH_TIME("BosonStarLevel::preCheckpointLevel");
    //Thomas Version for EMTENSOR
    /*fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(
                    MatterConstraints<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge(),
                  Density<ComplexScalarFieldWithPotential>(
                  complex_scalar_field, m_dx, m_p.G_Newton)),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);*/

     fillAllGhosts();
     Potential potential(m_p.potential_params);
     ComplexScalarFieldWithPotential complex_scalar_field(potential);
     BoxLoops::loop(make_compute_pack(
                     MatterConstraints<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge(),
                     EMTensor<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                     Interval(c_s11,c_s33))),
                     m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

}

// Things to do before outputting a plot file
void BosonStarLevel::prePlotLevel()
{
    CH_TIME("BosonStarLevel::prePlotLevel");
    //Thomas Version for EMTENSOR
    /*fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(
                    MatterConstraints<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge(),
                    Density<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton)),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);*/

     fillAllGhosts();
     Potential potential(m_p.potential_params);
     ComplexScalarFieldWithPotential complex_scalar_field(potential);
     BoxLoops::loop(make_compute_pack(
                     MatterConstraints<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge(),
                     EMTensor<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                     Interval(c_s11,c_s33))),
                     m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

}

// Things to do in RHS update, at each RK4 step
void BosonStarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
        a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ComplexScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    MatterCCZ4<ComplexScalarFieldWithPotential> my_ccz4_matter(
        complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Pi_Im + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_matter, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void BosonStarLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Things to do for analysis after each timestep and at the start
void BosonStarLevel::doAnalysis()
{
    CH_TIME("BosonStarLevel::specificPostTimeStep");
    bool first_step = (m_time == 0.0);
    if (m_p.activate_mass_extraction == 1)
    {
        // First compute the ADM Mass integrand values on the grid
        fillAllGhosts();
        auto weyl4_adm_compute_pack =
            make_compute_pack(Weyl4(m_p.extraction_params.extraction_center, m_dx),
                              ADMMass(m_p.L, m_dx));
        BoxLoops::loop(weyl4_adm_compute_pack, m_state_new, m_state_new,
                        EXCLUDE_GHOST_CELLS);

        // Do the extraction on the min extraction level
        {
          if (m_level == m_p.extraction_params.min_extraction_level())
            if (m_verbosity)
            {
                pout() << "BinaryBSLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }

            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }

        /* // Do the extraction on the min extraction level
        if (m_level == m_p.mass_extraction_params.min_extraction_level)
        {
            if (m_verbosity)
                pout() << "Extracting Mass." << endl;

            // Now refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            MassExtraction mass_extraction(m_p.mass_extraction_params, m_dt,
                                        m_time, first_step, m_restart_time);
            mass_extraction.execute_query(m_gr_amr.m_interpolator);
        } */
    }

    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(
                    MatterConstraints<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge()),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    if (m_level == 0)
    {
        if (m_p.calculate_constraint_violations)
        {
            // Write constraint violations to file
            ConstraintViolations constraint_violations(c_Ham,
                Interval(c_Mom1, c_Mom3), &m_gr_amr, m_p.coarsest_dx, m_dt,
                m_time, m_restart_time, "ConstraintViolations",
                first_step);
            constraint_violations.execute();
        }

        if (m_p.calculate_noether_charge)
        {
            // compute integrated volume weighted noether charge integral
            double noether_charge = m_gr_amr.compute_sum(c_N, m_dx);
            SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
                                            m_restart_time,
                                            SmallDataIO::APPEND,
                                            first_step);
            noether_charge_file.remove_duplicate_time_data();
            if (m_time == 0.)
            {
                noether_charge_file.write_header_line({"Noether Charge"});
            }
            noether_charge_file.write_time_data_line({noether_charge});
        }

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = m_gr_amr.compute_max(
                                Interval(c_mod_phi, c_mod_phi));
        SmallDataIO mod_phi_max_file("mod_phi_max", m_dt, m_time,
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});

        // Calculate the infinity-norm of all variables specified in params file
        // and output them
        if (m_p.num_vars_inf_norm > 0)
        {
            pout() << "Variable infinity norms:\n";
            for (int icomp : m_p.vars_inf_norm)
            {
                std::string var_name = UserVariables::variable_names[icomp];
                double var_norm = m_gr_amr.compute_norm(Interval(icomp, icomp),
                                    0., m_p.coarsest_dx);
                pout() << var_name << ": " << var_norm << "\t";
            }
            pout() << std::endl;
        }
    }


    if (m_p.gaussfit_params.do_star_tracking && m_level==m_p.gaussfit_params.AMR_level)
    {
        GaussianFitTracking gaussian_fit_tracking(m_p.gaussfit_params,m_dt,
                                        m_time,m_restart_time,first_step,m_p.L,m_level);

        gaussian_fit_tracking.do_star_tracking(m_gr_amr.m_interpolator);
        std::vector<double> dummy;
        gaussian_fit_tracking.get_BH_centres(dummy);
    }

    double S_phi_integral; // integral of angmomsource
    std::vector<double> S_phi_integrals(m_p.angmomflux_params.num_extraction_radii); // vector storing all integrals
    //if (m_p.do_flux_integration && m_level==m_p.angmomflux_params.extraction_level)
    if (m_p.do_flux_integration && m_level==m_p.angmomflux_params.min_extraction_level())
    {
        //std::cout << "Level : " << m_level << std::endl;
        // update stress tensor and mom flux components
        BoxLoops::loop(EMTensor_and_mom_flux<ComplexScalarFieldWithPotential>(
                      complex_scalar_field, m_dx, m_p.L, m_p.angmomflux_params.center,
                      c_Fphi_flux, c_Sphi_source, c_rho, Interval(c_s1,c_s3),
                      Interval(c_s11,c_s33)),  m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

        for (int i=m_p.angmomflux_params.num_extraction_radii-1; i>=0; i--)
        {
            // set angmomsource to zero outside of extraction radii
            BoxLoops::loop(SourceIntPreconditioner<ComplexScalarFieldWithPotential>(
                          complex_scalar_field, m_dx, m_p.L, m_p.angmomflux_params.center,
                          c_Sphi_source, m_p.angmomflux_params.extraction_radii[i]),
                          m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
            S_phi_integral = m_gr_amr.compute_sum(c_Sphi_source, m_dx);
            S_phi_integrals[i] = S_phi_integral;
            //std::cout << "t : " << m_time << ", i : " << i <<
            //          ", S_phi_integral : " << S_phi_integral << std::endl;



        }

        // save the Source integral to dat file
        std::vector<string> title_line(m_p.angmomflux_params.num_extraction_radii);
        string dummy_string;
        for (int j=0; j<m_p.angmomflux_params.num_extraction_radii; j++)
        {
            dummy_string = "r = " + to_string(m_p.angmomflux_params.extraction_radii[j]);
            title_line[j] = dummy_string;
        }

        SmallDataIO flux_file("AngMomSource", m_dt, m_time,
                                      m_restart_time,
                                      SmallDataIO::APPEND,
                                      first_step);

        if (m_time > 0) flux_file.remove_duplicate_time_data();

        if (m_time == 0.)
        {
            flux_file.write_header_line(title_line);
        }

        flux_file.write_time_data_line(S_phi_integrals);


        // Refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        // setup and do angmomflux integral
        AngMomFlux ang_mom_flux(m_p.angmomflux_params,m_time,m_dt,m_restart_time,first_step);
        ang_mom_flux.run(m_gr_amr.m_interpolator);
    }
}

// Specify if you want any plot files to be written, with which vars
void BosonStarLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = m_p.plot_vars;
}

void BosonStarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                   m_p.mass_extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi), current_state, tagging_criterion);
}
