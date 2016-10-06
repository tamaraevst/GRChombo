#if !defined(MATTERSFLEVEL_HPP_)
#error "This file should only be included through MatterSFLevel.hpp"
#endif

#ifndef MATTERSFLEVEL_IMPL_HPP_
#define MATTERSFLEVEL_IMPL_HPP_

#include "MatterSFLevel.hpp"
#include "FABDriver.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"
#include "ConstraintsMatter.hpp"
#include "CCZ4SFMatter.hpp"
#include "RelaxationChi.hpp"

//Initial data
#include "BubbleSF.hpp"

void MatterSFLevel::specificAdvance()
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) FABDriver<NanCheck>().execute(m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

void MatterSFLevel::initialData()
{
    CH_TIME("MatterSFLevel::initialData");
    if (m_verbosity) pout () << "MatterSFLevel::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc
    m_state_new.setVal(0.);

    //Initial conditions for scalar field - here a bubble
    FABDriver<BubbleSF>(m_p.sfm_params, m_dx).execute(m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());

}

void MatterSFLevel::preCheckpointLevel()
{
    fillAllGhosts();
    FABDriver<ConstraintsMatter>(m_dx).execute(m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void MatterSFLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{

    //Relaxation function for chi - this will eventually be done separately with hdf5 as input
    if (m_time < m_p.relaxtime) {

       //Calculate chi relaxation right hand side
       //Note that this assumes conformal chi and that the momentum constraint is trivially satisfied
       FABDriver<RelaxationChi>(m_dx, m_p.relaxspeed).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

       //No evolution in other variables, which are assumed to satisfy constraints per initial conditions
       a_rhs.setVal(0., Interval(c_h11,c_Mom3));
    } 

    //Else do normal CCZ4 evolution
    else {

    	FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

   	 	//Enforce positive chi and alpha
    	FABDriver<PositiveChiAndAlpha>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

    	//Calculate CCZ4 right hand side with SF matter
    	FABDriver<CCZ4SFMatter>(m_p.ccz4Params, m_dx, m_p.sigma).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

    	//We don't want undefined values floating around in the constraints
    	a_rhs.setVal(0., Interval(c_Ham,c_Mom3));

    }
}

void MatterSFLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);
}

#endif
