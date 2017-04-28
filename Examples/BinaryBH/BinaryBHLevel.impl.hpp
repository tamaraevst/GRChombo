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

//Initial data
#include "BinaryBH.hpp"

void BinaryBHLevel::specificAdvance()
{
    //Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(std::make_tuple(EnforceTfA(), PositiveChiAndAlpha()), m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    //if (m_p.nan_check) BoxLoops::loop<NanCheck>().execute(m_state_new, m_state_new, SKIP_GHOST_CELLS, no_simd_support());
#warning no_simd_support currently doesn't work! TODO: fix this!
}

void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity) pout () << "BinaryBHLevel::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc
#warning: TODO: Could write trivial compute class for this so that it can be included with other things
    m_state_new.setVal(0.);

    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);
    BoxLoops::loop(binary, m_state_new, m_state_new, FILL_GHOST_CELLS);
}

void BinaryBHLevel::preCheckpointLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void BinaryBHLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{
    //Enforce positive chi and alpha and trace free A
    BoxLoops::loop(std::make_tuple(EnforceTfA(), PositiveChiAndAlpha()), a_soln, a_soln, FILL_GHOST_CELLS);

    //Calculate CCZ4 right hand side
    BoxLoops::loop(CCZ4(m_p.ccz4Params, m_dx, m_p.sigma), a_soln, a_rhs, SKIP_GHOST_CELLS);

    //We don't want undefined values floating around in the constraints
    a_rhs.setVal(0., Interval(c_Ham,c_Mom3));
}

void BinaryBHLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free A_ij condition
    BoxLoops::loop(EnforceTfA(), a_soln, a_soln, FILL_GHOST_CELLS);
}

#endif
