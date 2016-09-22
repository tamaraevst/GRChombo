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
    CCZ4LevelFactory(ProfilingInfo * profilingInfo=nullptr);

    virtual
    AMRLevel* new_amrlevel() const;

    virtual
    ~CCZ4LevelFactory();

protected:
    Real m_dt_multiplier;
    SimulationParameters m_p;
    int m_tagBufferSize;
    int m_verbosity;
    ProfilingInfo* m_profilingInfo;
};


CCZ4LevelFactory::CCZ4LevelFactory (ProfilingInfo * a_profilingInfo):
    m_profilingInfo (a_profilingInfo)
{
    ParmParse pp;
    m_p.readParams(pp);
    pp.get("dt_multiplier", m_dt_multiplier);
    pp.get("tag_buffer_size", m_tagBufferSize);
    pp.get("verbosity", m_verbosity);
}
CCZ4LevelFactory::~CCZ4LevelFactory ()
{
}

// "virtual constructor"
AMRLevel*
CCZ4LevelFactory::new_amrlevel() const
{
    CCZ4Level* ccz4_level_ptr = new CCZ4Level (m_p, m_tagBufferSize, m_verbosity, m_profilingInfo);
    ccz4_level_ptr->initialDtMultiplier(m_dt_multiplier);
    return (static_cast <AMRLevel*> (ccz4_level_ptr));
}
#endif /* CCZ4LEVELFACTORY_HPP_ */
