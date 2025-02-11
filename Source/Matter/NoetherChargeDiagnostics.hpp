/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NOETHERCHARGEDIAGNOSTICS_HPP_
#define NOETHERCHARGEDIAGNOSTICS_HPP_

#include "ComplexScalarField.hpp"
#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
template <class deriv_t = FourthOrderDerivatives>
class NoetherChargeDiagnostics
{
protected:

    deriv_t m_deriv;

    // Need matter variables and chi
    template <class data_t> using ADMVars
                                = ADMConformalVars::VarsNoGauge<data_t>;
    template <class data_t> using MatterVars
                                = ComplexScalarField<>::Vars<data_t>;
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

public:

    NoetherChargeDiagnostics(const double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto vars = current_cell.template load_vars<Vars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto advec_csf =
            this->m_deriv.template advection<MatterVars>(current_cell, vars.shift);
        const auto advec =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        using namespace TensorAlgebra;

        // calculate Noether charge
        data_t N = pow(vars.chi, -1.5) * (matter_vars.phi_Im
            * matter_vars.Pi_Re - matter_vars.phi_Re * matter_vars.Pi_Im);

        data_t mod_phi = sqrt(matter_vars.phi_Re * matter_vars.phi_Re
                            + matter_vars.phi_Im * matter_vars.phi_Im);

        //Calculate time derivative phi^2
        data_t phi_Re_t = advec_csf.phi_Re - vars.lapse * matter_vars.Pi_Re;
        data_t phi_Im_t = advec_csf.phi_Im - vars.lapse * matter_vars.Pi_Im;
        data_t dt_A_sq = abs(2 * phi_Re_t * matter_vars.phi_Re + 2 * phi_Im_t * matter_vars.phi_Im);
        
        //Calculate time derivative of g_{tt}
        data_t lapse_t = 1.0 * advec.lapse - 2.0 * pow(vars.lapse, 1.0) * (vars.K - 2. * vars.Theta);
        Tensor<1, data_t> shift_t;
        FOR(i)
        {
            shift_t[i] = 0.75 * vars.B[i];
        } 

        data_t abs_gamma_tt;
        Tensor<2, data_t> h_t;

        data_t divshift = compute_trace(d1.shift);
        FOR(i, j)
        {
            h_t[i][j] = advec.h[i][j] - 2.0 * vars.lapse * vars.A[i][j] -
                      (2.0 / GR_SPACEDIM) * vars.h[i][j] * divshift;
            FOR(k)
            {
                h_t[i][j] +=
                vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
            }
        }

        data_t gamma_tt = -2.0 * vars.lapse * lapse_t; 
        FOR(i, j)
        {
            gamma_tt += shift_t[i] * vars.h[i][j] * vars.shift[j] + vars.shift[i] * (h_t[i][j] * shift_t[j] + vars.h[i][j] * shift_t[j]);
        }

        abs_gamma_tt = abs(gamma_tt);

        current_cell.store_vars(N, c_N);
        current_cell.store_vars(mod_phi, c_mod_phi);
        current_cell.store_vars(dt_A_sq, c_dt_mod_phi);
        current_cell.store_vars(abs_gamma_tt, c_gamma_tt);
    }
};

#endif /* NOETHERCHARGEDIAGNOSTICS_HPP_ */
