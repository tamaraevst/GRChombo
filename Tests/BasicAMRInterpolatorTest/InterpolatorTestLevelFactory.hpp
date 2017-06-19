#ifndef INTERPOLATORTESTLEVELFACTORY_HPP_
#define INTERPOLATORTESTLEVELFACTORY_HPP_

//General includes:
#include "AMRLevelFactory.H"
#include "ProfilingInfo.hpp"
#include "SimulationParameters.hpp"

//Problem specific includes:
#include "InterpolatorTestLevel.hpp"

class InterpolatorTestLevelFactory : public AMRLevelFactory
{
public:
    InterpolatorTestLevelFactory(SimulationParameters& a_sim_params);

    virtual
    AMRLevel* new_amrlevel() const;

    virtual
    ~InterpolatorTestLevelFactory();

protected:
    SimulationParameters m_p;
};


InterpolatorTestLevelFactory::InterpolatorTestLevelFactory (SimulationParameters& a_sim_params):
    m_p (a_sim_params)
{
}

InterpolatorTestLevelFactory::~InterpolatorTestLevelFactory ()
{
}

// "virtual constructor"
AMRLevel*
InterpolatorTestLevelFactory::new_amrlevel() const
{
    InterpolatorTestLevel* interpolator_test_level_ptr = new InterpolatorTestLevel (m_p, 0);
//    interpolator_test_level_ptr->initialDtMultiplier(m_p.dt_multiplier);
    return (static_cast <AMRLevel*> (interpolator_test_level_ptr));
}
#endif /* INTERPOLATORTESTLEVELFACTORY_HPP_ */
