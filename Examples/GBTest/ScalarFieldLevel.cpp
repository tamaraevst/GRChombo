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
#include "MatterOnly.hpp"
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"
#include "ChiTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "DefaultPotential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "ComputeModifiedScalars.hpp"
#include "GBAnalyticScalar.hpp"
#include "DebuggingTools.hpp"
#include "Coordinates.hpp"
#include <iostream>

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
    
    SetValue my_scalar_data(m_p.amplitude_scalar, Interval(c_phi, c_phi));

    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          my_scalar_data),
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
    
    BoxLoops::loop(make_compute_pack(
        Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3))),
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
    DefaultPotential potential;
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
    // Fixed grid - no pre-tagging
    // Pre tagging - fill ghost cells and calculate Ham terms
    //fillAllEvolutionGhosts();
    //BoxLoops::loop(make_compute_pack(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3), c_Ham_abs_terms, Interval(c_Moms_abs_terms, c_Moms_abs_terms))),
   //     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
   // We only use chi in the tagging criterion so only fill the ghosts for chi
   // fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
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

    fillAllGhosts();
    BoxLoops::loop(GBAnalyticScalar(m_dx, m_p.kerr_params.center, m_p.inner_r, m_p.outer_r),
                                            m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd());

    if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);

   bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' is
                        // called during setup at t=0 from Main
    
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

