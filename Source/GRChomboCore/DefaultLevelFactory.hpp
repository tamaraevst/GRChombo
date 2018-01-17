#ifndef DEFAULTLEVELFACTORY_HPP_
#define DEFAULTLEVELFACTORY_HPP_

// General includes:
#include "AMRLevelFactory.H"
#include "ProfilingInfo.hpp"
#include "SimulationParameters.hpp"

template <class level_t> class DefaultLevelFactory : public AMRLevelFactory
{
  public:
    DefaultLevelFactory(SimulationParameters &a_sim_params,
                        ProfilingInfo *a_profiling_info = nullptr)
        : m_p(a_sim_params), m_profiling_info(a_profiling_info)
    {
    }

    virtual AMRLevel *new_amrlevel() const
    {
        level_t *level_ptr = new level_t(m_p, m_p.verbosity, m_profiling_info);
        level_ptr->initialDtMultiplier(m_p.dt_multiplier);
        return (static_cast<AMRLevel *>(level_ptr));
    }

    virtual ~DefaultLevelFactory() {}

  protected:
    SimulationParameters m_p;
    ProfilingInfo *m_profiling_info;
};
#endif /* DEFAULTLEVELFACTORY_HPP_ */
