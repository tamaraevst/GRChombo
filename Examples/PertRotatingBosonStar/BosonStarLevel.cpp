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
// #include "IntegratedMovingPunctureGauge.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"
#include "NewConstraints.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"
#include "MovingBoxesRefinement.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "RotatingBosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

// For mass extraction
#include "ADMMass.hpp"
//#include "Density.hpp"
#include "EMTensor.hpp"
#include "MomFluxCalc.hpp"
#include "SourceIntPreconditioner.hpp"
#include "ADMMassExtraction.hpp"

// For GW extraction
#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"

// For Noether Charge calculation
#include "SmallDataIO.hpp"
#include "NoetherChargeDiagnostics.hpp"
#include "ModeDecomposition.hpp"

// For Ang Mom Integrating
#include "AngMomFlux.hpp"

// for chombo grid Functions
#include "AMRReductions.hpp"

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
    // if (m_verbosity)
        pout() << "BosonStarLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    RotatingBosonStar rotating_boson_star(m_p.rotating_bosonstar_params, m_dx);

    // read in BS profile 
    rotating_boson_star.compute_1d_rotating_solution();

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), rotating_boson_star),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(GammaCalculator(m_dx),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                   disable_simd());

    fillAllGhosts();
    // BoxLoops::loop(IntegratedMovingPunctureGauge(m_p.ccz4_params),
    //             m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a checkpoint file
void BosonStarLevel::preCheckpointLevel()
{
    CH_TIME("BosonStarLevel::preCheckpointLevel");

     fillAllGhosts();
     Potential potential(m_p.potential_params);
     ComplexScalarFieldWithPotential complex_scalar_field(potential);
     BoxLoops::loop(make_compute_pack(
                     MatterWeyl4<ComplexScalarFieldWithPotential>(
                     complex_scalar_field,m_p.extraction_params.extraction_center,
                     m_dx, m_p.formulation, m_p.G_Newton),
                     MatterConstraints<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                     Interval(c_Mom1, c_Mom3)), NoetherChargeDiagnostics<FourthOrderDerivatives>(m_dx),
                     EMTensor<ComplexScalarFieldWithPotential>(
                     complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                     Interval(c_s11,c_s33))),
                     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

}

// Things to do before outputting a plot file
void BosonStarLevel::prePlotLevel()
{
    CH_TIME("BosonStarLevel::prePlotLevel");

      fillAllGhosts();
      Potential potential(m_p.potential_params);
      ComplexScalarFieldWithPotential complex_scalar_field(potential);
      BoxLoops::loop(make_compute_pack(
                      MatterWeyl4<ComplexScalarFieldWithPotential>(
                      complex_scalar_field,m_p.extraction_params.extraction_center,
                      m_dx, m_p.formulation, m_p.G_Newton),
                      MatterConstraints<ComplexScalarFieldWithPotential>(
                      complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                      Interval(c_Mom1, c_Mom3)),  NoetherChargeDiagnostics<FourthOrderDerivatives>(m_dx),
                      EMTensor<ComplexScalarFieldWithPotential>(
                      complex_scalar_field, m_dx, c_rho, Interval(c_s1,c_s3),
                      Interval(c_s11,c_s33))),
                      m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

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
    MatterCCZ4RHS<ComplexScalarFieldWithPotential> my_ccz4_matter(
        complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    // MatterCCZ4RHS<ComplexScalarFieldWithPotential, IntegratedMovingPunctureGauge, FourthOrderDerivatives> my_ccz4_matter(
    //     complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
    //     m_p.G_Newton);
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

    // First compute the ADM Mass integrand values on the grid
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    auto weyl4_adm_compute_pack =
               make_compute_pack(MatterWeyl4<ComplexScalarFieldWithPotential>(
               complex_scalar_field,m_p.extraction_params.extraction_center,
               m_dx, m_p.formulation, m_p.G_Newton), ADMMass(m_p.center, m_dx));
    BoxLoops::loop(weyl4_adm_compute_pack, m_state_new, m_state_diagnostics,
                        EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(MatterConstraints<ComplexScalarFieldWithPotential>(
                        complex_scalar_field, m_dx, m_p.G_Newton, c_Ham,
                        Interval(c_Mom1, c_Mom3)), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    
    if (m_p.activate_weyl_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        CH_TIME("BosonStarLevel::doAnalysis::Weyl4&ADMMass");
        
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
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
    }

    if (m_p.activate_mass_extraction == 1 &&
        m_level == m_p.mass_extraction_params.min_extraction_level())
    {
        if (m_verbosity)
        {
            pout() << "BinaryBSLevel::specificPostTimeStep:"
                      " Extracting mass." << endl;
        }

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        ADMMassExtraction mass_extraction(m_p.mass_extraction_params, m_dt,
                                    m_time, first_step,
                                    m_restart_time);
        mass_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    // noether charge, max mod phi, min chi, constraint violations
    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherChargeDiagnostics<FourthOrderDerivatives>(m_dx), m_state_new, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
    }
    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
        if (m_p.calculate_noether_charge)
        {
            //noether charge should be calculated pre-check and pre plot
            //so automatically here

            // compute integrated volume weighted noether charge integral

            double noether_charge = amr_reductions.sum(c_N);
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
        
        double dtA_00_re = amr_reductions.sum(c_dt_mod_phi);
        double dtA_01_re = amr_reductions.sum(c_mode_dtA_01_re);
        double dtA_01_im = amr_reductions.sum(c_mode_dtA_01_im);
        double dtA_02_re = amr_reductions.sum(c_mode_dtA_02_re);
        double dtA_02_im = amr_reductions.sum(c_mode_dtA_02_im);
        double dtA_03_re = amr_reductions.sum(c_mode_dtA_03_re);
        double dtA_03_im = amr_reductions.sum(c_mode_dtA_03_im);
        double dtA_04_re = amr_reductions.sum(c_mode_dtA_04_re);
        double dtA_04_im = amr_reductions.sum(c_mode_dtA_04_im);
        double dtA_05_re = amr_reductions.sum(c_mode_dtA_05_re);
        double dtA_05_im = amr_reductions.sum(c_mode_dtA_05_im);
        double dtA_06_re = amr_reductions.sum(c_mode_dtA_06_re);
        double dtA_06_im = amr_reductions.sum(c_mode_dtA_06_im);
        double dtA_07_re = amr_reductions.sum(c_mode_dtA_07_re);
        double dtA_07_im = amr_reductions.sum(c_mode_dtA_07_im);
        double dtA_08_re = amr_reductions.sum(c_mode_dtA_08_re);
        double dtA_08_im = amr_reductions.sum(c_mode_dtA_08_im);
        double dtA_09_re = amr_reductions.sum(c_mode_dtA_09_re);
        double dtA_09_im = amr_reductions.sum(c_mode_dtA_09_im);
        double dtA_10_re = amr_reductions.sum(c_mode_dtA_10_re);
        double dtA_10_im = amr_reductions.sum(c_mode_dtA_10_im);
        double dtA_11_re = amr_reductions.sum(c_mode_dtA_11_re);
        double dtA_11_im = amr_reductions.sum(c_mode_dtA_11_im);
        double dtA_12_re = amr_reductions.sum(c_mode_dtA_12_re);
        double dtA_12_im = amr_reductions.sum(c_mode_dtA_12_im);
        double dtA_13_re = amr_reductions.sum(c_mode_dtA_13_re);
        double dtA_13_im = amr_reductions.sum(c_mode_dtA_13_im);
        double dtA_14_re = amr_reductions.sum(c_mode_dtA_14_re);
        double dtA_14_im = amr_reductions.sum(c_mode_dtA_14_im);
        double dtA_15_re = amr_reductions.sum(c_mode_dtA_15_re);
        double dtA_15_im = amr_reductions.sum(c_mode_dtA_15_im);
        double dtA_16_re = amr_reductions.sum(c_mode_dtA_16_re);
        double dtA_16_im = amr_reductions.sum(c_mode_dtA_16_im);
        double dtA_17_re = amr_reductions.sum(c_mode_dtA_17_re);
        double dtA_17_im = amr_reductions.sum(c_mode_dtA_17_im);
        double dtA_18_re = amr_reductions.sum(c_mode_dtA_18_re);
        double dtA_18_im = amr_reductions.sum(c_mode_dtA_18_im);
        double dtA_19_re = amr_reductions.sum(c_mode_dtA_19_re);
        double dtA_19_im = amr_reductions.sum(c_mode_dtA_19_im);
        double dtA_20_re = amr_reductions.sum(c_mode_dtA_20_re);
        double dtA_20_im = amr_reductions.sum(c_mode_dtA_20_im);

        SmallDataIO dtA_modes_file("dtA_modes", m_dt, m_time,
            m_restart_time,
            SmallDataIO::APPEND,
            first_step);
        dtA_modes_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            dtA_modes_file.write_header_line({"dtA00", "dtA01 Re", "dtA01 Im", "dtA02 Re", "dtA02 Im", "dtA03 Re", "dtA03 Im", "dtA04 Re", "dtA04 Im", "dtA05 Re", "dtA05 Im", "dtA06 Re", "dtA06 Im", "dtA07 Re", "dtA07 Im", "dtA08 Re", "dtA08 Im", "dtA09 Re", "dtA09 Im", "dtA10 Re", "dtA10 Im", "dtA11 Re", "dtA11 Im", "dtA12 Re", "dtA12 Im", "dtA13 Re", "dtA13 Im", "dtA14 Re", "dtA14 Im", "dtA15 Re", "dtA15 Im", "dtA16 Re", "dtA16 Im", "dtA17 Re", "dtA17 Im", "dtA18 Re", "dtA18 Im", "dtA19 Re", "dtA19 Im", "dtA20 Re", "dtA20 Im"});
        }
        dtA_modes_file.write_time_data_line({dtA_00_re, dtA_01_re, dtA_01_im, dtA_02_re, dtA_02_im, dtA_03_re, dtA_03_im, dtA_04_re, dtA_04_im, dtA_05_re, dtA_05_im, dtA_06_re, dtA_06_im, dtA_07_re, dtA_07_im, dtA_08_re, dtA_08_im, dtA_09_re, dtA_09_im, dtA_10_re, dtA_10_im, dtA_11_re, dtA_11_im, dtA_12_re, dtA_12_im, dtA_13_re, dtA_13_im, dtA_14_re, dtA_14_im, dtA_15_re, dtA_15_im, dtA_16_re, dtA_16_im, dtA_17_re, dtA_17_im, dtA_18_re, dtA_18_im, dtA_19_re, dtA_19_im, dtA_20_re, dtA_20_im});

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = amr_reductions.max(c_mod_phi);
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


        // Compute the min of chi and write it to a file
        AMRReductions<VariableType::evolution> amr_reductions_ev(m_gr_amr);
        double min_chi = amr_reductions_ev.min(c_chi);
        SmallDataIO min_chi_file("min_chi", m_dt, m_time,
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        min_chi_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            min_chi_file.write_header_line({"min chi"});
        }
        min_chi_file.write_time_data_line({min_chi});
        
	// Compute the min of lapse and write it to a file
        double min_lapse = amr_reductions_ev.min(c_lapse);
        SmallDataIO min_lapse_file("min_lapse", m_dt, m_time,
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        min_lapse_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            min_lapse_file.write_header_line({"min lapse"});
        }
        min_lapse_file.write_time_data_line({min_lapse});

        // constraeints calculated pre check and pre plot so done here already

        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 2, true);
        double L1_Ham = amr_reductions.norm(c_Ham, 1, true);
        double L1_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 1, true);
        SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                     m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom", "L^1_Ham", "L^1_Mom",});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom, L1_Ham, L1_Mom});

    // Compute the maximum of mod_phi and write it to a file
    double dt_phi_max = amr_reductions.max(c_dt_mod_phi);
    double dt_gamma_max = amr_reductions.max(c_gamma_tt);
    SmallDataIO dt_max_file("dt_max", m_dt, m_time,
                             m_restart_time,
                             SmallDataIO::APPEND,
                             first_step);
    dt_max_file.remove_duplicate_time_data();
    if (m_time == 0.)
    {
        dt_max_file.write_header_line({"max(dt phi)", "max(dt gtt)"});
    }
    dt_max_file.write_time_data_line({dt_phi_max, dt_gamma_max});
    }
}

void BosonStarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                    m_p.extraction_params, m_p.regrid_threshold_phi,
                    m_p.regrid_threshold_chi, m_p.activate_extraction), current_state, tagging_criterion);

//      BoxLoops::loop(MovingBoxesRefinement(
//                           m_dx, m_level, m_p.tag_puncture_max_level,
//                           m_p.center, m_p.puncture_radius, m_p.puncture_mass, m_p.tag_buffer),
//                      current_state, tagging_criterion);	

}
