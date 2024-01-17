/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STAMR_HPP_
#define STAMR_HPP_

#include "GRAMR.hpp"
#include "StarTracker.hpp"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for star spacetimes
 */
class STAMR : public GRAMR
{
  public:
    StarTracker m_star_tracker;

    STAMR() {}

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator) override
    {
        GRAMR::set_interpolator(a_interpolator);
        m_star_tracker.set_interpolator(a_interpolator);
    }
};

#endif /* STAMR_HPP_ */