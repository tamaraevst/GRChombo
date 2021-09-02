/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GBSCALARANALYTIC_HPP_
#define GBSCALARANALYTIC_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"

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

        data_t x = coords.x;
        data_t y = coords.y;
        data_t z = coords.z;

        const data_t r = coords.get_radius();

        double M = 1.0;
        double beta = 1.0;

        data_t xx = (1.0 / M) * r * pow(1 + M / (2.0 * r), 2.0);

        data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

        current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
};

#endif /* GBSCALARANALYTIC_HPP_ */
