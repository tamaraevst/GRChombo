#ifndef DEFAULTLEVELFACTORY_HPP_
#define DEFAULTLEVELFACTORY_HPP_

// General includes:
#include "AMRLevelFactory.H"
#include "SimulationParameters.hpp"

template <class level_t> class DefaultLevelFactory : public AMRLevelFactory
{
  public:
    DefaultLevelFactory(SimulationParameters &a_sim_params) : m_p(a_sim_params)
    {
    }

    virtual AMRLevel *new_amrlevel() const
    {
        level_t *level_ptr = new level_t(m_p, m_p.verbosity);
        level_ptr->initialDtMultiplier(m_p.dt_multiplier);
        return (static_cast<AMRLevel *>(level_ptr));
    }

    virtual ~DefaultLevelFactory() {}

  protected:
    SimulationParameters m_p;
};
#endif /* DEFAULTLEVELFACTORY_HPP_ */
