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

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "BosonStar.hpp"
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

#include "ComputeDiagnostics.hpp"
#include "MetricxxExtraction.hpp"
#include "MetricxyExtraction.hpp"
#include "MetricxzExtraction.hpp"
#include "MetricyyExtraction.hpp"
#include "MetricyzExtraction.hpp"
#include "MetriczzExtraction.hpp"

#include "MetricxxrExtraction.hpp"
#include "MetricxyrExtraction.hpp"
#include "MetricxzrExtraction.hpp"
#include "MetricyyrExtraction.hpp"
#include "MetricyzrExtraction.hpp"
#include "MetriczzrExtraction.hpp"

#include "MetricxxtExtraction.hpp"
#include "MetricxytExtraction.hpp"
#include "MetricxztExtraction.hpp"
#include "MetricyytExtraction.hpp"
#include "MetricyztExtraction.hpp"
#include "MetriczztExtraction.hpp"

#include "ShiftxExtraction.hpp"
#include "ShiftyExtraction.hpp"
#include "ShiftzExtraction.hpp"

#include "ShiftrxExtraction.hpp"
#include "ShiftryExtraction.hpp"
#include "ShiftrzExtraction.hpp"

#include "ShifttxExtraction.hpp"
#include "ShifttyExtraction.hpp"
#include "ShifttzExtraction.hpp"

#include "LapseExtraction.hpp"
#include "LapsetExtraction.hpp"
#include "LapserExtraction.hpp"

// For Noether Charge calculation
#include "SmallDataIO.hpp"
#include "NoetherCharge.hpp"

#include "ComputeWeightFunction.hpp"

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
    if (m_verbosity)
        pout() << "BosonStarLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    BosonStar boson_star(m_p.bosonstar_params, m_p.bosonstar2_params, m_p.potential_params,
                         m_p.G_Newton, m_dx, m_p.identical, m_verbosity);


    // the max radius the code might need to calculate out to is L*sqrt(3)
    boson_star.compute_1d_solution(4.*m_p.L);

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(GammaCalculator(m_dx),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                   disable_simd());

    BoxLoops::loop(ComputeWeightFunction(m_p.bosonstar_params, m_p.bosonstar2_params, m_dx), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd());

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
                     Interval(c_Mom1, c_Mom3)), NoetherCharge()),
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
                      Interval(c_Mom1, c_Mom3)), NoetherCharge()),
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
    BoxLoops::loop(ComputeDiagnostics<FourthOrderDerivatives>(m_p.ccz4_params.lapse_advec_coeff, m_p.ccz4_params.lapse_coeff, m_p.ccz4_params.lapse_power, m_p.ccz4_params.shift_advec_coeff, m_p.ccz4_params.shift_Gamma_coeff, m_dx, m_p.center), m_state_new, m_state_diagnostics,
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

            MetricxxExtraction gxx_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxx_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxyExtraction gxy_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxy_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxzExtraction gxz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxz_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyyExtraction gyy_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyy_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyzExtraction gyz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyz_extraction.execute_query(m_gr_amr.m_interpolator);

            MetriczzExtraction gzz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gzz_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxxrExtraction gxxr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxxr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxxtExtraction gxxt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxxt_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxyrExtraction gxyr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxyr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxytExtraction gxyt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxyt_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxzrExtraction gxzr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxzr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricxztExtraction gxzt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gxzt_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyyrExtraction gyyr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyyr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyytExtraction gyyt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyyt_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyzrExtraction gyzr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyzr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetricyztExtraction gyzt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gyzt_extraction.execute_query(m_gr_amr.m_interpolator);

            MetriczzrExtraction gzzr_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gzzr_extraction.execute_query(m_gr_amr.m_interpolator);

            MetriczztExtraction gzzt_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gzzt_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftxExtraction shiftx_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftx_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftyExtraction shifty_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shifty_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftzExtraction shiftz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftz_extraction.execute_query(m_gr_amr.m_interpolator);

            ShifttxExtraction shifttx_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shifttx_extraction.execute_query(m_gr_amr.m_interpolator);

            ShifttyExtraction shiftty_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftty_extraction.execute_query(m_gr_amr.m_interpolator);

            ShifttzExtraction shifttz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shifttz_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftrxExtraction shiftrx_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftrx_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftryExtraction shiftry_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftry_extraction.execute_query(m_gr_amr.m_interpolator);

            ShiftrzExtraction shiftrz_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            shiftrz_extraction.execute_query(m_gr_amr.m_interpolator);

            LapseExtraction lapse_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            lapse_extraction.execute_query(m_gr_amr.m_interpolator);

            LapserExtraction lapser_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            lapser_extraction.execute_query(m_gr_amr.m_interpolator);

            LapsetExtraction lapset_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            lapset_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }


    // noether charge, max mod phi, min chi, constraint violations
    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherCharge(), m_state_new, m_state_diagnostics,
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
        double min_chi = amr_reductions.min(c_chi);
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

        // constraints calculated pre check and pre plot so done here already

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
    }
}

void BosonStarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
   BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                   m_p.extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi, m_p.activate_extraction), current_state, tagging_criterion);

}