#if !defined(KERRBHLEVEL_HPP_)
#error "This file should only be included through KerrBHLevel.hpp"
#endif

#ifndef KERRBHLEVEL_IMPL_HPP_
#define KERRBHLEVEL_IMPL_HPP_

#include "KerrBHLevel.hpp"
#include "BoxLoops.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"
#include "Constraints.hpp"
#include "CCZ4.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"

//Initial data
#include "KerrBH.hpp"
#include "GammaCalculator.hpp"

void KerrBHLevel::specificAdvance()
{
    //Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS, no_simd_support());
}

void KerrBHLevel::initialData()
{
    CH_TIME("KerrBHLevel::initialData");
    if (m_verbosity) pout () << "KerrBHLevel::initialData " << m_level << endl;

    //First set everything to zero (to avoid undefinded values on constraints) then calculate initial data
    //Get the Kerr solution in the variables, then calculate the \tilde\Gamma^i numerically as these 
    //are non zero and not calculated in the Kerr ICs
    BoxLoops::loop(make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx)), 
                                                      m_state_new, m_state_new, FILL_GHOST_CELLS);

    fillAllGhosts(); 
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void KerrBHLevel::preCheckpointLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void KerrBHLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{
    //Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), a_soln, a_soln, FILL_GHOST_CELLS);


    //Calculate CCZ4 right hand side and set constraints to zero to avoid undefined values
    BoxLoops::loop(make_compute_pack(CCZ4(m_p.ccz4Params, m_dx, m_p.sigma), SetValue(0, Interval(c_Ham, c_Mom3)) ),
                   a_soln, a_rhs, SKIP_GHOST_CELLS);

    //For now this seems necessary
    a_rhs.setVal(0.0, Interval(c_Ham, c_Mom3));

}

void KerrBHLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free A_ij condition
    BoxLoops::loop(EnforceTfA(), a_soln, a_soln, FILL_GHOST_CELLS);
}

// Specify which variables to write at plot intervals
void KerrBHLevel::specificWritePlotHeader(std::vector<int>& plot_states) const
{
    //Specify the variables we want to output as plot
    plot_states = {c_chi, c_K, c_lapse, c_shift1};    
}

#endif
