/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GBSCALARANALYTIC_HPP_
#define GBSCALARANALYTIC_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "CoordinateTransformations.hpp"
#include "UserVariables.hpp"
#include "CoordinateTransformations.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
#include <cmath>
#include "CCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"
#include "FourthOrderDerivatives.hpp"

template <class gauge_t = MovingPunctureGauge, class deriv_t = FourthOrderDerivatives>
class GBScalarAnalytic
{
  public:

    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;
    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    GBScalarAnalytic(std::array<double, CH_SPACEDIM> a_center, double a_dx, double a_mass, double a_beta_amplitude) : m_center(a_center), m_dx(a_dx), m_mass(a_mass), m_beta_amplitude(a_beta_amplitude)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        //load vars
        const auto vars_is = current_cell.template load_vars<CCZ4Vars>();

        //where am I?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        
        //rename some parameters for convenience
        double M = m_mass;
        double beta = m_beta_amplitude;
       
        // Transform from Schwarzschild to isotropic and make conformal 
        data_t r =  vars_is.chi * coords.get_radius();
        data_t xx = pow((1.0 + M / (2.0 * r)), 2.0) * r / M;
        data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));
  
        current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_mass;
    const double m_beta_amplitude;

};

#endif /* GBSCALARANALYTIC_HPP_ */
