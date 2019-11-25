/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

#include "AMR.H"
#include "AMRInterpolator.hpp"
#include "Lagrange.hpp"
#include <chrono>
#include <ratio>

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy (such as the AMRInterpolator)
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */
class GRAMR : public AMR
{
  private:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();

    // the info for the puncture tracks
    int m_num_punctures;
    std::vector<std::array<double, CH_SPACEDIM>> m_puncture_coords;

  public:
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    GRAMR()
    {
        m_interpolator = nullptr;
        m_num_punctures = 2; // default to 2 for now
        m_puncture_coords.resize(m_num_punctures);
    }

    // function to set punctures
    void set_puncture_coords(
        std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords)
    {
        m_puncture_coords = a_puncture_coords;
    }

    // function to get punctures
    std::vector<std::array<double, CH_SPACEDIM>> get_puncture_coords()
    {
        return m_puncture_coords;
    }

    auto get_walltime()
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_time);

        return duration.count();
    }

    // Called after AMR object set up
    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }
};

#endif /* GRAMR_HPP_ */
