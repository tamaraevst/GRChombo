/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#include "PunctureTracker.hpp"

#if defined(USE_AHFINDER) && !defined(DISABLE_AHFINDER)
#include "AHInterpolation.hpp"
#include "AHSphericalCoords.hpp"
#include "ApparentHorizon.hpp"
#endif

// Chombo namespace
#include "UsingNamespace.H"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */
class BHAMR : public GRAMR
{
  public:
    PunctureTracker m_puncture_tracker;

    BHAMR() {}

   ~BHAMR()
    {
#if defined(USE_AHFINDER) && !defined(DISABLE_AHFINDER)
        // destroy horizon pointers and finalize PETSc
        for (auto &ah : m_apparent_horizons)
            delete ah;
        AHFinder::finalize();
#endif
    }

#if defined(USE_AHFINDER) && !defined(DISABLE_AHFINDER)
    std::vector<ApparentHorizon<AHInterpolation<AHSphericalCoords>> *>
        m_apparent_horizons;
#endif

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator) override
    {
        GRAMR::set_interpolator(a_interpolator);
        m_puncture_tracker.set_interpolator(a_interpolator);
    }
};

#endif /* BHAMR_HPP_ */
