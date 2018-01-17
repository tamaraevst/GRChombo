#ifndef KERRBHLEVEL_HPP_
#define KERRBHLEVEL_HPP_

#include "GRAMRLevel.hpp"

class KerrBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<KerrBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    virtual void specificAdvance();

    // initialize data
    virtual void initialData();

    virtual void preCheckpointLevel();

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs, Real a_dt);

    // Specify which variables to write at plot intervals
    virtual void specificWritePlotHeader(std::vector<int> &plot_states) const;

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#include "KerrBHLevel.impl.hpp"

#endif /* KERRBHLEVEL_HPP_ */
