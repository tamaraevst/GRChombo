#if !defined(BINARYBHLEVEL_HPP_)
#error "This file should only be included through BinaryBHLevel.hpp"
#endif

#ifndef BINARYBHLEVEL_IMPL_HPP_
#define BINARYBHLEVEL_IMPL_HPP_

#include "BinaryBHLevel.hpp"
#include "BoxLoops.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"
#include "Constraints.hpp"
#include "CCZ4.hpp"
#include "ComputeClassPack.hpp"
#include "SetValue.hpp"

//Initial data
#include "BinaryBH.hpp"

void BinaryBHLevel::specificAdvance()
{
    //Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS, no_simd_support());
}

void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity) pout () << "BinaryBHLevel::initialData " << m_level << endl;

    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx); //Set up the compute class for the BinaryBH initial data
    //First set everything to zero (to avoid undefinded values on constraints) then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new, m_state_new, FILL_GHOST_CELLS);
}

void BinaryBHLevel::preCheckpointLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void BinaryBHLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{
    //Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), a_soln, a_soln, FILL_GHOST_CELLS);

    //Calculate CCZ4 right hand side and set constraints to zero to avoid undefined values
    BoxLoops::loop(make_compute_pack( CCZ4(m_p.ccz4Params, m_dx, m_p.sigma), SetValue(0, Interval(c_Ham, c_Mom3)) ),
                   a_soln, a_rhs, SKIP_GHOST_CELLS);
}

void BinaryBHLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free A_ij condition
    BoxLoops::loop(EnforceTfA(), a_soln, a_soln, FILL_GHOST_CELLS);
}

#endif
