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
#include "ADMConformalVars.hpp"
#include "simd.hpp"
#include <cmath>
#include "DebuggingTools.hpp"
#include "KerrBH.hpp"

class GBScalarAnalytic
{
  public:
    /// BSSN variables
    template <class data_t> using MetricVars = BSSNVars::VarsWithGauge<data_t>;
    GBScalarAnalytic(std::array<double, CH_SPACEDIM> a_center, double a_dx, double a_mass, double a_beta_amplitude) : m_center(a_center), m_dx(a_dx), m_mass(a_mass), m_beta_amplitude(a_beta_amplitude)
    {
    }
    // using params_t = typename KerrBH::params_t;

    // template <class data_t>
    // using MetricVars = typename KerrBH::template Vars<data_t>;

    // GBScalarAnalytic(params_t a_params, double a_dx, double a_beta_amplitude) : KerrBH(a_params, a_dx), m_beta_amplitude(a_beta_amplitude)
    // {
    // }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        const auto vars_is = current_cell.template load_vars<MetricVars>();

        // Tensor<2, data_t> spherical_g;
        // MetricVars<data_t> vars_is;
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // compute_isotropic_metric(spherical_g, coords);

        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        
        double M = m_mass;
        double beta = m_beta_amplitude;
        DEBUG_OUT(vars_is.chi);
        // vars_is.h = CoordinateTransformations::spherical_to_cartesian_LL(spherical_g, x, y, z);

        // using namespace TensorAlgebra;
        // // // Convert to BSSN vars
        // data_t deth = compute_determinant(vars_is.h);
        // auto h_UU = compute_inverse_sym(vars_is.h);
        // vars_is.chi = pow(deth, -1. / 3.);

        //For PG coordinates
        // data_t r = (1.0 / vars.chi) * coords.get_radius();
        // data_t xx = r / M;
        // data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

        // For isotropic coordinates
        data_t r =  coords.get_radius();
        data_t xx = sqrt(vars_is.chi) * pow((1 + M / (2.0 * r)), 2) * r / M;
        data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));
        DEBUG_OUT(xx);

        DEBUG_OUT(phi_analytic);
        current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_mass;
    const double m_beta_amplitude;

};

#endif /* GBSCALARANALYTIC_HPP_ */
