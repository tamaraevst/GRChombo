#ifndef MATTERSFLEVEL_HPP_
#define MATTERSFLEVEL_HPP_

#include "GRAMRLevel.hpp"

class MatterSFLevel : public GRAMRLevel
{
    friend class MatterSFLevelFactory;
    //Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    virtual
    void specificAdvance();

    // initialize data
    virtual
    void initialData();

    virtual
    void preCheckpointLevel();

    virtual
    void specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time);

    virtual
    void specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt);
};

#include "MatterSFLevel.impl.hpp"

#endif /* MATTERSFLEVEL_HPP_ */
