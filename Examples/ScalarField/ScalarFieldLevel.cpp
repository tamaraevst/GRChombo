/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

#include "FilesystemTools.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "DefaultPotential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "ComputeModifiedScalars.hpp"
#include "GBScalarAnalytic.hpp"

// For post processing
#include "SmallDataIO.hpp"
#include "AMRReductions.hpp"


// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    DefaultPotential potential();
    ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);

    BoxLoops::loop(make_compute_pack(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom))),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    DefaultPotential potential();
    ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
    
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(
        FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
        current_state, tagging_criterion);
}

//Output norms of Gauss-Bonnet and Chern-Simons into file 
void ScalarFieldLevel::specificPostTimeStep()
{
    CH_TIME("ScalarFieldLevel::specificPostTimeStep");
    DefaultPotential potential();
    ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);
    
    if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);

    // bool first_step = (m_time == 0.);
    bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.calculate_scalar_norm)
    {
        fillAllGhosts();
        BoxLoops::loop(ComputeModifiedScalars(m_p.center, m_dx,
                     m_p.gamma_amplitude, 
                     m_p.beta_amplitude),
                     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        if (m_level ==0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double CS_norm = amr_reductions.norm(c_chernsimons, 1, true); // L1 norm of Chern Simons
            double GB_norm = amr_reductions.norm(c_gaussbonnet, 1, true); // L1 norm of Gauss Bonnet

            if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);
            SmallDataIO scalars_file(m_p.data_path + "modified_scalars_l1norm",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            scalars_file.remove_duplicate_time_data();
            if (first_step)
                {
                    scalars_file.write_header_line({"norm_ChernSimons", "norm_GaussBonnet"});
                }
            scalars_file.write_time_data_line({CS_norm, GB_norm});
        }
    }

      if (m_p.compare_gb_analytic)
    {
        fillAllGhosts();
        BoxLoops::loop(GBScalarAnalytic(m_p.center, m_dx), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_red_diag(m_gr_amr);
            AMRReductions<VariableType::evolution> amr_red_ev(m_gr_amr);

            //output norms
            double DiffWithAnalytic = amr_red_diag.norm(c_phinumerical - c_phianalytic, 1, true);
            // double NormPhi = amr_red_ev.norm(c_phi, 1, true);
            SmallDataIO norm_phi_file(m_p.data_path + "normdiff_phi_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_phi_file.remove_duplicate_time_data();
            if (first_step)
                {
                    // norm_phi_file.write_header_line({"Phi Norm", "Analytic Phi Norm"});
                    norm_phi_file.write_header_line({"Phi Error Norm"});
                }
            // norm_phi_file.write_time_data_line({NormPhi, NormAnalytic});
            norm_phi_file.write_time_data_line({DiffWithAnalytic});
        }
    }
}

