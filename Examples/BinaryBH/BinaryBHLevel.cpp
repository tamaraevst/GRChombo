
/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "BinaryBH.hpp"
#include "ScalarField.hpp"
#include "BoxLoops.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "BinaryPunctureTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "NewConstraints.hpp"
#include "MatterCCZ4RHS.hpp"
#include "FixedGridsTaggingCriterion.hpp"
#include "ChiAndPhiTaggingCriterion.hpp"
#include "TraceARemoval.hpp"
// #include "InitialScalarData.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"


#include "DefaultPotential.hpp"
#include "PunctureTracker.hpp"

#include "SetValue.hpp"
#include "SmallDataIO.hpp"

#include "AMRReductions.hpp"
#include "ComputeModifiedScalars.hpp"

#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"
#include "PhiExtraction.hpp"


// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBHLevel::initialData " << m_level << endl;
#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // set the value of phi - constant over the grid
    SetValue my_scalar_data(m_p.amplitude_scalar, Interval(c_phi, c_phi));

    // Can't use simd with this initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.),two_punctures_initial_data, my_scalar_data), m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);
    // set the value of phi - constant over the grid
    SetValue my_scalar_data(m_p.amplitude_scalar, Interval(c_phi, c_phi));

    // First set everything to zero (to avoid undefinded values)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary, my_scalar_data),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
#endif
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    DefaultPotential potential;
    ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);
    MatterCCZ4RHS<ScalarFieldWithPotential> my_ccz4_matter(
        scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// void BinaryBHLevel::preTagCells()
// {
//     // We only use chi in the tagging criterion so only fill the ghosts for chi
//     fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
// }

// specify the cells to tag
void BinaryBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                            const FArrayBox &current_state)
{
    std::array<double, 2> puncture_masses;
    std::vector<std::array<double, CH_SPACEDIM>> punctures;

    if (m_p.track_punctures)
    {
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
        punctures = m_tp_amr.m_puncture_tracker.get_puncture_coords();
    }
    BoxLoops::loop(BinaryPunctureTaggingCriterion<FourthOrderDerivatives>(
                       m_dx, m_level, m_p.tag_horizons_max_levels,
                       m_p.tag_punctures_max_levels, m_p.extraction_params,
                       punctures, m_p.activate_extraction, m_p.track_punctures,
                       puncture_masses, m_p.bh_tagging_buffers,
                       m_p.puncture_tag_min_separation),
                   current_state, tagging_criterion);

    // set tagging criterion for scalar field
    BoxLoops::loop(
            ChiAndPhiTaggingCriterion(
            m_dx,  m_p.regrid_threshold_chi, m_p.regrid_threshold_phi),
            current_state, tagging_criterion, disable_simd());
    
}

void BinaryBHLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHLevel::specificPostTimeStep");
    if (m_verbosity)
        pout() << "BinaryBHLevel::specificPostTimestep " << m_level << endl;

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main
    AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.activate_extraction_phi == 1)
    {   
        int min_level = m_p.extraction_params_phi.min_extraction_level();
        bool calculate_phi = at_level_timestep_multiple(min_level);
        if (calculate_phi)
        {
            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("PhiExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::evolution, Interval(c_phi, c_Pi_Re),
                    min_level);
                PhiExtraction my_extraction(m_p.extraction_params_phi, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {   
        if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);
            
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)), m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                         m_restart_time, SmallDataIO::APPEND,
                                         first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }

    if (m_p.calculate_scalar_norm)
    {   
        if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);
            
        fillAllGhosts();
        BoxLoops::loop(ComputeModifiedScalars(m_p.center, m_dx,
                     m_p.gamma_amplitude, 
                     m_p.beta_amplitude),
                     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        if (m_level ==0)
        {
            double CS_norm = amr_reductions.norm(c_chernsimons, 1, true); // L1 norm of Chern Simons
            double GB_norm = amr_reductions.norm(c_gaussbonnet, 1, true); // L1 norm of Gauss Bonnet

            SmallDataIO scalars_file(m_p.data_path + "modified_scalars_l1norm",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            scalars_file.remove_duplicate_time_data();
            if (first_step)
                {
                    scalars_file.write_header_line({"Norm Chern Simons", "Norm Gauss Bonnet"});
                }
            scalars_file.write_time_data_line({CS_norm, GB_norm});
        }
    }

    if (m_p.compute_all_norms)
    {   
        if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);
            
        fillAllGhosts();

        if (m_level == 0)
        {   
            AMRReductions<VariableType::evolution> amr_red_ev(m_gr_amr);
            //output norms
            double NormNumericPhi = amr_red_ev.norm(c_phi, 1, true);

            double NormChi = amr_red_ev.norm(c_chi, 1, true);
            double NormK = amr_red_ev.norm(c_K, 1, true);

            double Normh11 = amr_red_ev.norm(c_h11, 1, true);
            double Normh13 = amr_red_ev.norm(c_h13, 1, true);
            double Normh12 = amr_red_ev.norm(c_h12, 1, true);
            double Normh22 = amr_red_ev.norm(c_h22, 1, true);
            double Normh23 = amr_red_ev.norm(c_h23, 1, true);
            double Normh33 = amr_red_ev.norm(c_h33, 1, true);

            double NormA11 = amr_red_ev.norm(c_A11, 1, true);
            double NormA12 = amr_red_ev.norm(c_A11, 1, true);
            double NormA13 = amr_red_ev.norm(c_A13, 1, true);
            double NormA22 = amr_red_ev.norm(c_A22, 1, true);
            double NormA23 = amr_red_ev.norm(c_A23, 1, true);
            double NormA33 = amr_red_ev.norm(c_A33, 1, true);
        

            SmallDataIO norm_phi_file(m_p.data_path + "norm_phi_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_phi_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_phi_file.write_header_line({"Phi Numeric Norm"});
                }
            norm_phi_file.write_time_data_line({NormNumericPhi});

            SmallDataIO norm_chiK_file(m_p.data_path + "norm_chiK_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_chiK_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_chiK_file.write_header_line({"Chi", "K"});
                }
            norm_chiK_file.write_time_data_line({NormChi, NormK});

            SmallDataIO norm_A_file(m_p.data_path + "norm_A_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_A_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_A_file.write_header_line({"A11", "A12", "A13", "A22", "A23", "A33"});
                }
            norm_A_file.write_time_data_line({NormA11, NormA12, NormA13, NormA22, NormA23, NormA33});

            SmallDataIO norm_h_file(m_p.data_path + "norm_h_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_h_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_h_file.write_header_line({"h11", "h12", "h13", "h22", "h23", "h33"});
                }
            norm_h_file.write_time_data_line({Normh11, Normh12, Normh13, Normh22, Normh23, Normh33});

        }
    }
}

#ifdef CH_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars
void BinaryBHLevel::prePlotLevel()
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::prePlotLevel" << m_level << endl;
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3))),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
    
}
#endif /* CH_USE_HDF5 */
