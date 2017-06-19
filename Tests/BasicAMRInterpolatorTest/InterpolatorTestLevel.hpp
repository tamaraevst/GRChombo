#ifndef INTERPOLATORTESTLEVEL_HPP_
#define INTERPOLATORTESTLEVEL_HPP_

#include "GRAMRLevel.hpp"

class InterpolatorTestLevel : public GRAMRLevel
{
    friend class InterpolatorTestLevelFactory;
    //Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual
    void initialData()
    {
        m_state_new.setVal(42.);
    }

    virtual
    void specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time) {}
};

#endif /* INTERPOLATORTESTLEVEL_HPP_ */
