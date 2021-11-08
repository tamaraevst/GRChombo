/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "DebuggingTools.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionDiagnostics
{
  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;
    const double m_inner_r;
    const double m_outer_r;

  public:
    ExcisionDiagnostics(const double a_dx,
                        const std::array<double, CH_SPACEDIM> a_center,
                        background_t a_background, const double a_inner_r,
                        const double a_outer_r)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center),
          m_background(a_background), m_inner_r(a_inner_r), m_outer_r(a_outer_r)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        double horizon_distance = m_background.excise(current_cell);
        if (coords.get_radius() < m_inner_r || coords.get_radius() > m_outer_r)
        {
	    if (coords.get_radius() > m_outer_r)
	    { 
		double r = coords.get_radius();
		DEBUG_OUT(r);
		double xx = pow((1.0 + 1.0 / (2.0 * r)), 2.0) * r;
                double phi = (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));
	        DEBUG_OUT(phi);
	    }

	    if (coords.get_radius() < m_inner_r)
            {
                double r_in = coords.get_radius();
                DEBUG_OUT(r_in);
                double xx_in = pow((1.0 + 1.0 / (2.0 * r_in)), 2.0) * r_in;
                double phi_in = (1.0 / xx_in + 1.0 / (xx_in * xx_in) + (4.0 / 3.0) * 1.0 / (xx_in * xx_in * xx_in));
                DEBUG_OUT(phi_in);
            }
            current_cell.store_vars(0.0, c_xMom);
            current_cell.store_vars(0.0, c_phianalytic);
        } // else do nothing
    }
};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
