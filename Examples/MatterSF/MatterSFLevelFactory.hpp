#ifndef MATTERSFLEVELFACTORY_HPP_
#define MATTERSFLEVELFACTORY_HPP_

//General includes:
#include "AMRLevelFactory.H"
#include "ProfilingInfo.hpp"
#include "SimulationParameters.hpp"

//Problem specific includes:
#include "MatterSFLevel.hpp"

class MatterSFLevelFactory : public AMRLevelFactory
{
public:
    MatterSFLevelFactory(SimulationParameters& a_sim_params, ProfilingInfo * profilingInfo=nullptr);

    virtual
    AMRLevel* new_amrlevel() const;

    virtual
    ~MatterSFLevelFactory();

protected:
    Real m_dt_multiplier;
    SimulationParameters m_p;
    ProfilingInfo* m_profilingInfo;
};


MatterSFLevelFactory::MatterSFLevelFactory (SimulationParameters& a_sim_params, ProfilingInfo * a_profilingInfo):
    m_p (a_sim_params), m_profilingInfo (a_profilingInfo)
{
}

MatterSFLevelFactory::~MatterSFLevelFactory ()
{
}

// "virtual constructor"
AMRLevel*
MatterSFLevelFactory::new_amrlevel() const
{
    MatterSFLevel* matter_sf_level_ptr = new MatterSFLevel (m_p, m_p.verbosity, m_profilingInfo);
    matter_sf_level_ptr->initialDtMultiplier(m_p.dt_multiplier);
    return (static_cast <AMRLevel*> (matter_sf_level_ptr));
}
#endif /* MATTERSFLEVELFACTORY_HPP_ */
