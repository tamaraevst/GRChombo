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

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "ComplexPhiAndKTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "BosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

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
    BosonStar boson_star(m_p.bosonstar_params, m_p.potential_params, m_dx);

    // the max radius the code might need to calculate out to is 2*sqrt(3)*L
    // 3.5 is an upper bound to 2*sqrt(3)
    const double max_radius{m_p.L*3.5};

    boson_star.compute_1d_solution(max_radius);

    /*
    //For debugging
    pout() << setprecision(12);
    pout() << setw(16) << "R" << "\t" << setw(16) << "chi" << "\t" <<
    setw(16) << "lapse" << "\t" << setw(16) << "phi" "\n";
    double r;
    for(int i = 1; i < 2 * static_cast<int>(m_p.L/m_dx); ++i)
    {
        r = i * m_dx;
        pout() << setw(16) << r << "\t" << setw(16) <<
        boson_star.m_1d_sol.m_chi(r) << "\t" << setw(16) <<
        boson_star.m_1d_sol.m_lapse(r) << "\t" << setw(16) <<
        boson_star.m_1d_sol.m_phi(r) << "\n";
    }
    */

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), boson_star),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());
}

// Things to do before outputting a checkpoint file
void BosonStarLevel::preCheckpointLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ComplexScalarFieldWithPotential>(
                       complex_scalar_field, m_dx, m_p.G_Newton),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void BosonStarLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(MatterConstraints<ComplexScalarFieldWithPotential>(
                       complex_scalar_field, m_dx, m_p.G_Newton),
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
    SetValue set_constraints_zero(0.0, Interval(c_Ham, c_Mom3));
    auto compute_pack2 =
        make_compute_pack(my_ccz4_matter, set_constraints_zero);
    BoxLoops::loop(compute_pack2, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

}

// Things to do at ODE update, after soln + rhs
void BosonStarLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Specify if you want any plot files to be written, with which vars
void BosonStarLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = {c_phi_Re, c_phi_Im, c_chi, c_Ham};
}

void BosonStarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ComplexPhiAndKTaggingCriterion(m_dx, m_p.regrid_threshold_phi,
                                           m_p.regrid_threshold_K),
                   current_state, tagging_criterion);
}
