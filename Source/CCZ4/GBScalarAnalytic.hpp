/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GBSCALARANALYTIC_HPP_
#define GBSCALARANALYTIC_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"
#include "CoordinateTransformations.hpp"
#include <cmath>
class GBScalarAnalytic
{
  public:
    GBScalarAnalytic(std::array<double, CH_SPACEDIM> &a_center, double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        Tensor<2, data_t> spherical_g;

        compute_PL_metric(spherical_g, coords);

        data_t area_element = CoordinateTransformations::area_element_sphere(spherical_g);

        data_t ar_radius = area_element / (4 * M_PI);

        data_t areal_radius = pow(ar_radius, 0.5);
        
        double M = 1.0;

        double beta = 1.0;

        data_t xx = areal_radius / M;

        data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

        current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;

    template <class data_t> void compute_PL_metric(Tensor<2, data_t> &spherical_g, 
                                        const Coordinates<data_t> coords) const
    {
        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;

        // the radius, subject to a floor
        data_t r = coords.get_radius();
        data_t r2 = r * r;

        // calculate useful position quantities
        data_t rho2 = simd_max(x * x + y * y, 1e-12);
        data_t rho = sqrt(rho2);
        data_t sin_theta = rho / r;
        data_t sin_theta2 = sin_theta * sin_theta;

        // Metric in PL coordinates
        spherical_g[0][0] = 1; 
        spherical_g[0][1] = 0;
        spherical_g[0][2] = 0;

        spherical_g[1][0] = 0;
        spherical_g[1][1] = r2;
        spherical_g[1][2] = 0;

        spherical_g[2][0] = 0;
        spherical_g[2][1] = 0;
        spherical_g[2][2] = r2 * sin_theta2;

    }

};

#endif /* GBSCALARANALYTIC_HPP_ */
