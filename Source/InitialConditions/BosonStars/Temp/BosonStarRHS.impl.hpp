/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARRHS_HPP_)
#error "This file should only be included through BosonStarRHS.hpp"
#endif

#ifndef BOSONSTARRHS_IMPL_HPP_
#define BOSONSTARRHS_IMPL_HPP_

BosonStarRHS::BosonStarRHS(Potential::params_t a_params_potential,
    double a_G_Newton = 1.0)
    : m_rescaled_phi4_coeff(a_params_potential.phi4_coeff
        / (4.0 * M_PI * a_G_Newton * a_params_potential.scalar_mass *
        a_params_potential.scalar_mass)), m_G_Newton(a_G_Newton) {}

template <typename initial_state_t>
void BosonStarRHS::operator() (const initial_state_t &a_vars,
initial_state_t &rhs_out, const double &a_rho)
{
    //introduce references to make equations easier to read
    auto &alpha = a_vars[0]; //alpha = (1/2)*log(g_tt) - frequency/m
    auto &beta = a_vars[1]; //beta = (1/2)*log(g_rr)
    auto &psi = a_vars[2]; //psi = sqrt(4*pi*G)*|phi|
    auto &Psi = a_vars[3]; //Psi = d(psi)/d(rho)

    auto &alpha_rhs = rhs_out[0];
    auto &beta_rhs = rhs_out[1];
    auto &psi_rhs = rhs_out[2];
    auto &Psi_rhs = rhs_out[3];

    //the RHS contains removable singularities so need to handle the case
    //rho = 0 separately
    if( a_rho < std::numeric_limits<double>::epsilon() )
    {
        alpha_rhs = 0.0;
        beta_rhs = 0.0;
        psi_rhs = 0.0;
        Psi_rhs = 0.0;
    }
    else
    {
        //first calculate some quantities which are used more than once
        auto exp_minus2alpha = std::exp(-2.0 * alpha);
        auto exp_plus2beta = std::exp(2.0 * beta);
        auto psi2 = psi * psi;
        auto Psi2 = Psi * Psi;
        auto rho_inv = 1.0 / a_rho;

        alpha_rhs = 0.5 * (a_rho * (exp_plus2beta * psi2 * (exp_minus2alpha - 1.0
            - 0.5 * m_rescaled_phi4_coeff * psi2) + Psi2)
            + rho_inv * (exp_plus2beta - 1.0));
        beta_rhs = 0.5 * (a_rho * (exp_plus2beta * psi2 * (exp_minus2alpha + 1.0
            + 0.5 * m_rescaled_phi4_coeff * psi2) + Psi2)
            - rho_inv * (exp_plus2beta -1.0));
        psi_rhs = Psi;
        Psi_rhs = (exp_plus2beta * a_rho * psi2 * (1.0
            + 0.5 * m_rescaled_phi4_coeff * psi2) - rho_inv * (exp_plus2beta +
            1.0)) * Psi + exp_plus2beta * (1.0 + m_rescaled_phi4_coeff * psi2 -
            exp_minus2alpha) * psi;
    }
}

#endif /* BOSONSTARRHS_IMPL_HPP_ */
