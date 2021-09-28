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
#include "ComplexScalarPotential.hpp"
#include "ExcisionDiagnostics.hpp"
#include "ExcisionEvolution.hpp"
#include "FixedBGScalarField.hpp"
#include "FixedBGEvolution.hpp"
#include "InitialConditions.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "DebuggingTools.hpp"

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
    InitialConditions set_phi(m_p.field_amplitude, m_p.center,
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
        scalar_field, boosted_bh, m_p.sigma, m_dx, m_p.center);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

    // excise within horizon, no simd
    BoxLoops::loop(
        ExcisionEvolution<ScalarFieldWithPotential, IsotropicKerrFixedBG>(
            m_dx, m_p.center, boosted_bh),
        a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
}

void ScalarFieldLevel::specificPostTimeStep()
{
    // At any level, but after the coarsest timestep
    double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;
    const double remainder = fmod(m_time, coarsest_dt);
    if (min(abs(remainder), abs(remainder - coarsest_dt)) < 1.0e-8)
    {
        // calculate the density of the PF, but excise the BH region completely
        fillAllGhosts();
        DefaultPotential potential;
        ScalarFieldWithPotential scalar_field(potential);
        IsotropicKerrFixedBG boosted_bh(m_p.bg_params, m_dx);
        // excise within horizon
        BoxLoops::loop(
            ExcisionDiagnostics<ScalarFieldWithPotential, IsotropicKerrFixedBG>(
                m_dx, m_p.center, boosted_bh, m_p.inner_r, m_p.outer_r),
            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS,
            disable_simd());

        AMRReductions<VariableType::diagnostic> amr_red_diag(m_gr_amr);
        if (m_level == 0)
        {
            bool first_step = (m_time == m_dt);
            double NormAnalyticPhi = amr_red_diag.norm(c_phianalytic, 1); //we are not normalising over volume here so we need to divide by appropriate
                                                                            //volume post-processing
            SmallDataIO norm_phi_analytic_file("norm_phi_analytic_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_phi_analytic_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_phi_analytic_file.write_header_line({"Phi Analytic Norm"});
                }
            norm_phi_analytic_file.write_time_data_line({NormAnalyticPhi});
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
