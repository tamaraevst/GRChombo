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
// #include "BSSNVars.hpp"
#include "simd.hpp"
#include <cmath>

class GBScalarAnalytic
{
  public:
    /// BSSN variables
    template <class data_t> using Vars = BSSNVars::VarsWithGauge<data_t>;
    
    GBScalarAnalytic(std::array<double, CH_SPACEDIM> &a_center, double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // const auto vars = current_cell.template load_vars<Vars>();

        Tensor<2, data_t> spherical_g;
        Vars<data_t> vars;
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        compute_isotropic_metric(spherical_g, coords);

        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        
        double M = 1.0;
        double beta = 1.0;
        
        vars.h = CoordinateTransformations::spherical_to_cartesian_LL(spherical_g, x, y, z);

        using namespace TensorAlgebra;
        // Convert to BSSN vars
        data_t deth = compute_determinant(vars.h);
        auto h_UU = compute_inverse_sym(vars.h);
        vars.chi = pow(deth, -1. / 3.);

        //For PG coordinates
        // data_t r = (1.0 / vars.chi) * coords.get_radius();
        // data_t xx = r / M;
        // data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

        // For isotropic coordinates
        data_t r_conformal = vars.chi * coords.get_radius();
        data_t xx = pow((1 + M / (2.0 * r_conformal)), 2) * r_conformal / M;
        data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

        current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;

    template <class data_t> void compute_isotropic_metric(Tensor<2, data_t> &spherical_g, 
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

        double M = 1.0;

        // Metric in isotropic coordinates
        spherical_g[0][0] = pow((1 + M / (2.0 * r)), 4); 
        spherical_g[0][1] = 0;
        spherical_g[0][2] = 0;

        spherical_g[1][0] = 0;
        spherical_g[1][1] = r2 * pow((1 + M / (2.0 * r)), 2);
        spherical_g[1][2] = 0;

        spherical_g[2][0] = 0;
        spherical_g[2][1] = 0;
        spherical_g[2][2] = r2 * sin_theta2 * pow((1 + M / (2.0 * r)), 2);

    }

};

#endif /* GBSCALARANALYTIC_HPP_ */
