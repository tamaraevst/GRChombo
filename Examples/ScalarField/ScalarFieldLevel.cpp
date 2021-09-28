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

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "DefaultPotential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "ComputeModifiedScalars.hpp"

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
    DefaultPotential potential;
    ScalarFieldWithPotential scalar_field(potential, m_p.gamma_amplitude, m_p.beta_amplitude);

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
    if (m_p.matter_only == true)
    {
        if (m_p.max_spatial_derivative_order == 4)
        {   
            MatterOnly<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
            BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

        }
        else if (m_p.max_spatial_derivative_order == 6)
        {   
            MatterOnly<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
            BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
        }
    }
    if (m_p.matter_only == false)
    {
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

    AMRReductions<VariableType::diagnostic> amr_red_diag(m_gr_amr);
    AMRReductions<VariableType::evolution> amr_red_ev(m_gr_amr);

    if (!FilesystemTools::directory_exists(m_p.data_path))
            FilesystemTools::mkdir_recursive(m_p.data_path);

   bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' is
                        // called during setup at t=0 from Main
    
    //Calculates L2 norms of Hamiltonian and Momentum constraints
    if (m_p.calculate_constraint_norms)
    {
        if (m_level == 0)
        {
            double L2_Ham = amr_red_diag.norm(c_Ham);
            double L2_Mom = amr_red_diag.norm(Interval(c_Mom1, c_Mom3));
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

    //Calculates L1 norms of GB and CS scalars
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

    //Calculates L1 norm of numeric \phi
    if (m_p.norm_numeric_phi)
    {   
        if (m_level == 0)
        {
            double NormNumericPhi = amr_red_ev.norm(c_phi, 1, true);

            SmallDataIO norm_phi_file(m_p.data_path + "norm_numeric_phi_values",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            norm_phi_file.remove_duplicate_time_data();
            if (first_step)
                {
                    norm_phi_file.write_header_line({"Phi Numeric Norm"});
                }
            norm_phi_file.write_time_data_line({NormNumericPhi});
        }
    }

    //Calculates L1 norms of geometric quantities, useful sometimes for postprocessing/figuring out what's going on
    if (m_p.geometric_norms)
    {
        if(m_level == 0)
        {   
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

