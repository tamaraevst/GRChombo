#ifndef CCZ4LEVEL_HPP_
#define CCZ4LEVEL_HPP_

#include "GRAMRLevel.hpp"

class BinaryBHLevel : public GRAMRLevel
{
    friend class BinaryBHLevelFactory;
    //Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    virtual
    void specificAdvance();

    virtual
    void specificPostTimeStep();

    // initialize data
    virtual
    void initialData();

    virtual
    void fillBdyGhosts();

    virtual
    void preCheckpointLevel();

    virtual
    void specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time);

    virtual
    void specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt);
};

#include "BinaryBHLevel.impl.hpp"

#endif /* CCZ4LEVEL_HPP_ */
