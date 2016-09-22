#ifndef CCZ4LEVEL_HPP_
#define CCZ4LEVEL_HPP_

#include "GRAMRLevel.hpp"

class CCZ4Level : public GRAMRLevel
{
    friend class CCZ4LevelFactory;
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

#endif /* CCZ4LEVEL_HPP_ */
