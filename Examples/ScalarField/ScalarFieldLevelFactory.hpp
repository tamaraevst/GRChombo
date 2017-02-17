#ifndef SCALARFIELDLEVELFACTORY_HPP_
#define SCALARFIELDLEVELFACTORY_HPP_

//General includes:
#include "AMRLevelFactory.H"
#include "ProfilingInfo.hpp"
#include "SimulationParameters.hpp"

//Problem specific includes:
#include "ScalarFieldLevel.hpp"

class ScalarFieldLevelFactory : public AMRLevelFactory
{
public:
    ScalarFieldLevelFactory(SimulationParameters& a_sim_params, ProfilingInfo * profilingInfo=nullptr);

    virtual
    AMRLevel* new_amrlevel() const;

    virtual
    ~ScalarFieldLevelFactory();

protected:
    Real m_dt_multiplier;
    SimulationParameters m_p;
    ProfilingInfo* m_profilingInfo;
};


ScalarFieldLevelFactory::ScalarFieldLevelFactory (SimulationParameters& a_sim_params, ProfilingInfo * a_profilingInfo):
    m_p (a_sim_params), m_profilingInfo (a_profilingInfo)
{
}

ScalarFieldLevelFactory::~ScalarFieldLevelFactory ()
{
}

// "virtual constructor"
AMRLevel*
ScalarFieldLevelFactory::new_amrlevel() const
{
    ScalarFieldLevel* matter_sf_level_ptr = new ScalarFieldLevel (m_p, m_p.verbosity, m_profilingInfo);
    matter_sf_level_ptr->initialDtMultiplier(m_p.dt_multiplier);
    return (static_cast <AMRLevel*> (matter_sf_level_ptr));
}
#endif /* SCALARFIELDLEVELFACTORY_HPP_ */
