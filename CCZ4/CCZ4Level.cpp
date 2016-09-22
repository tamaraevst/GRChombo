#include "CCZ4Level.hpp"
#include "FABDriver.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"
#include "Constraints.hpp"

//Initial data
#include "BinaryBH.hpp"

void CCZ4Level::specificAdvance()
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(m_state_new, m_state_new, true);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(m_state_new, m_state_new, true);

    //Check for nan's
    if (m_p.nan_check) FABDriver<NanCheck>().execute(m_state_new, m_state_new, false, disable_simd());
}

void CCZ4Level::specificPostTimeStep()
{
}

void CCZ4Level::initialData()
{
    CH_TIME("CCZ4Level::initialData");
    if (m_verbosity) pout () << "CCZ4Level::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc
    m_state_new.setVal(0.);

    FABDriver<BinaryBH>(m_p.bh1_params, m_p.bh2_params, m_dx).execute(m_state_new, m_state_new, true, disable_simd());
}

void CCZ4Level::fillBdyGhosts()
{
}

void CCZ4Level::preCheckpointLevel()
{
    fillAllGhosts();
    FABDriver<Constraints>(m_dx).execute(m_state_new, m_state_new, false);
}

void CCZ4Level::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{
    FABDriver<EnforceTfA>().execute(a_soln, a_soln, true);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(a_soln, a_soln, true);

    //Calculate CCZ4 right hand side
    FABDriver<CCZ4>(m_p.ccz4Params, m_dx, m_p.sigma).execute(a_soln, a_rhs, false);

    //We don't want undefined values floating around in the constraints
    a_rhs.setVal(0., Interval(c_Ham,c_Mom3));
}

void CCZ4Level::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(a_soln, a_soln, true);
}
