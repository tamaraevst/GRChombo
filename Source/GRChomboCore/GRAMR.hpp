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
#include <sys/time.h>

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy (such as the AMRInterpolator)
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */
class GRAMR : public AMR
{
  public:
    AMRInterpolator<Lagrange<4>> *m_interpolator; //!< The interpolator pointer

    GRAMR() { m_interpolator = nullptr; } // constructor
    ~GRAMR() { delete m_interpolator; }   // destructor

    auto get_walltime()
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_time);

        return duration.count();
    }

    void set_interpolator(const int dx_scalar, const int &a_verbosity = 0)
    {
        // Setup the AMRInterpolator
        std::array<double, CH_SPACEDIM> origin, dx;
        dx.fill(dx_scalar);
        origin.fill(dx_scalar / 2);
        AMRInterpolator<Lagrange<4>> interpolator(*this, origin, dx,
                                                  a_verbosity);
        m_interpolator = &interpolator;
    }

  private:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();
};

#endif /* GRAMR_HPP_ */
