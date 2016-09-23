#ifndef CCZ4LEVELFACTORY_HPP_
#define CCZ4LEVELFACTORY_HPP_

//General includes:
#include "AMRLevelFactory.H"
#include "ProfilingInfo.hpp"
#include "SimulationParameters.hpp"

//Problem specific includes:
#include "CCZ4Level.hpp"

class CCZ4LevelFactory : public AMRLevelFactory
{
public:
    CCZ4LevelFactory(SimulationParameters& a_sim_params, ProfilingInfo * profilingInfo=nullptr);

    virtual
    AMRLevel* new_amrlevel() const;

    virtual
    ~CCZ4LevelFactory();

protected:
    Real m_dt_multiplier;
    SimulationParameters m_p;
    ProfilingInfo* m_profilingInfo;
};


CCZ4LevelFactory::CCZ4LevelFactory (SimulationParameters& a_sim_params, ProfilingInfo * a_profilingInfo):
    m_p (a_sim_params), m_profilingInfo (a_profilingInfo)
{
}

CCZ4LevelFactory::~CCZ4LevelFactory ()
{
}

// "virtual constructor"
AMRLevel*
CCZ4LevelFactory::new_amrlevel() const
{
    CCZ4Level* ccz4_level_ptr = new CCZ4Level (m_p, m_p.verbosity, m_profilingInfo);
    ccz4_level_ptr->initialDtMultiplier(m_p.dt_multiplier);
    return (static_cast <AMRLevel*> (ccz4_level_ptr));
}
#endif /* CCZ4LEVELFACTORY_HPP_ */
