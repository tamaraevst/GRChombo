/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GBSCALARANALYTIC_HPP_
#define GBSCALARANALYTIC_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "ADMConformalVars.hpp"
#include "CoordinateTransformations.hpp"
#include "UserVariables.hpp"
#include "CoordinateTransformations.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
#include <cmath>
#include "DebuggingTools.hpp"

class GBScalarAnalytic
{
  public:

  template <class data_t>
  using Vars = ADMConformalVars::VarsWithGauge<data_t>;

    GBScalarAnalytic(std::array<double, CH_SPACEDIM> a_center, double a_dx, double a_mass, double a_beta_amplitude) : m_center(a_center), m_dx(a_dx), m_mass(a_mass), m_beta_amplitude(a_beta_amplitude)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // //load vars
        // const auto vars_is = current_cell.template load_vars<CCZ4Vars>();

      // The cartesian variables 
      Vars<data_t> vars;
      Tensor<2, data_t> spherical_g;

      //where am I?
      Coordinates<data_t> coords(current_cell, m_dx, m_center);
      // Compute the components in spherical coords as per 1401.1548
      compute_isotropic_metric(spherical_g, coords);

      // work out where we are on the grid
      data_t x = coords.x;
      double y = coords.y;
      double z = coords.z;

      using namespace CoordinateTransformations;
      // Convert spherical components to cartesian components using coordinate
      // transforms
      vars.h = spherical_to_cartesian_LL(spherical_g, x, y, z);
  
      using namespace TensorAlgebra;
      // Convert to BSSN vars
      data_t deth = compute_determinant(vars.h);
      auto h_UU = compute_inverse_sym(vars.h);
      vars.chi = pow(deth, -1. / 3.);
      DEBUG_OUT(vars.chi);

      //rename for convenience
      double M = m_mass;
      double beta = m_beta_amplitude;
       
      // Transform from Schwarzschild to isotropic and make conformal 
      data_t r =  sqrt(vars.chi) * coords.get_radius();
      data_t xx = pow((1.0 + M / (2.0 * r)), 2.0) * r / M;
      data_t phi_analytic = (2.0 * beta) / (M * M) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

      DEBUG_OUT(r);
      DEBUG_OUT(xx);
      DEBUG_OUT(phi_analytic);
  
      current_cell.store_vars(phi_analytic, c_phianalytic);
    }

   protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_mass;
    const double m_beta_amplitude;


    template <class data_t> void compute_isotropic_metric(Tensor<2, data_t> &spherical_g, 
                                        const Coordinates<data_t> coords) const
    {
        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        // the radius, subject to a floor
        data_t r = coords.get_radius();
        data_t r2 = r * r;

        // calculate useful position quantities
        data_t rho2 = simd_max(x * x + y * y, 1e-12);
        data_t rho = sqrt(rho2);
        data_t sin_theta = rho / r;
        data_t cos_theta = z / r;
        data_t sin_theta2 = sin_theta * sin_theta;

        double M = m_mass;

         // calculate useful metric quantities
        double r_plus = M + sqrt(M * M);
        double r_minus = M - sqrt(M * M);

        // The Boyer-Lindquist coordinate
        data_t r_BL = r * pow(1.0 + 0.25 * r_plus / r, 2.0);

        // Other useful quantities per 1001.4077
        data_t Sigma = r_BL * r_BL;
        // In the paper this is just 'A', but not to be confused with A_ij
        data_t AA = pow(r_BL * r_BL, 2.0);
        // The rr component of the conformal spatial matric
        data_t gamma_rr =
        Sigma * pow(r + 0.25 * r_plus, 2.0) / (r * r2 * (r_BL - r_minus));

        // Metric in semi isotropic Kerr-Schild coordinates, r, theta (t or th), phi
        // (p)
        FOR(i, j) { spherical_g[i][j] = 0.0; }
        spherical_g[0][0] = gamma_rr;                // gamma_rr
        spherical_g[1][1] = Sigma;                   // gamma_tt
        spherical_g[2][2] = AA / Sigma * sin_theta2; // gamma_p

    }


};

#endif /* GBSCALARANALYTIC_HPP_ */
