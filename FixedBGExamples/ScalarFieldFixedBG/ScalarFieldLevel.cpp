/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "DefaultPotential.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGScalarField.hpp"
#include "FixedBGEvolution.hpp"
#include "InitialConditions.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "ComputeModifiedScalars.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for fields
    SetValue set_zero(0.0);
    IsotropicKerrFixedBG boosted_bh(m_p.bg_params, m_dx); // just calculates chi
    InitialConditions set_phi(m_p.field_amplitude,
                              m_p.center,
                              m_p.bg_params, m_dx);
    auto compute_pack = make_compute_pack(set_zero, boosted_bh);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);
    BoxLoops::loop(set_phi, m_state_new, m_state_new, FILL_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, IsotropicKerrFixedBG>(
            m_dx, m_p.center, boosted_bh),
        m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a plot file
void ScalarFieldLevel::prePlotLevel() {}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    DefaultPotential potential;
    ScalarFieldWithPotential scalar_field(potential);
    IsotropicKerrFixedBG boosted_bh(m_p.bg_params, m_dx);
    FixedBGEvolution<ScalarFieldWithPotential, IsotropicKerrFixedBG> my_matter(
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center, m_p.beta_amplitude, m_p.gamma_amplitude);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, IsotropicKerrFixedBG>(
            m_dx, m_p.center, boosted_bh),
        a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::specificPostTimeStep()
{
    CH_TIME("ScalarFieldLevel::specificPostTimeStep");

    bool first_step = (m_time == m_dt); // if not called in Main

    AMRReductions<VariableType::diagnostic> amr_red_diag(m_gr_amr);
    AMRReductions<VariableType::evolution> amr_red_ev(m_gr_amr);
    // DefaultPotential potential;
    // ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);
    
    if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);

//    bool first_step =
//         (m_time == 0.); // this form is used when 'specificPostTimeStep' is
                        // called during setup at t=0 from Main

    if (m_p.calculate_scalar_norm)
    {   
        fillAllGhosts();
        BoxLoops::loop(ComputeModifiedScalars(m_p.center, m_dx,
                     m_p.gamma_amplitude, 
                     m_p.beta_amplitude),
                     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        if (m_level ==0)
        {
            double CS_norm = amr_red_diag.norm(c_chernsimons, 1, true); // L1 norm of Chern Simons
            double GB_norm = amr_red_diag.norm(c_gaussbonnet, 1, true); // L1 norm of Gauss Bonnet
            // double GB_norm_2 = amr_reductions.norm(c_gaussbonnet_2, 1, true); // L1 norm of Gauss Bonnet

            if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);
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

    if (m_p.compare_gb_analytic)
    {
        fillAllGhosts();

        if (m_level == 0)
        {
            //output norms
            double NormNumericPhi = amr_red_ev.norm(c_phi, 1, true);
            double NormAnalyticPhi = amr_red_diag.norm(c_phianalytic, 1, true);

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
                    norm_phi_file.write_header_line({"Phi Analytic Norm", "Phi Numeric Norm"});
                }
            norm_phi_file.write_time_data_line({NormAnalyticPhi, NormNumericPhi});

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

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.regrid_length,
                                              m_p.center),
                   current_state, tagging_criterion, disable_simd());
}
