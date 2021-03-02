/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MOMFLUXCALC_HPP)
#error "This file should only be included through MomFluxCalc.hpp"
#endif

#ifndef MOMFLUXCALC_IMPL_HPP
#define MOMFLUXCALC_IMPL_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "simd.hpp"


template <class matter_t>
EMTensor_and_mom_flux<matter_t>::EMTensor_and_mom_flux(const matter_t &a_matter,
                             const double dx, const double a_L,
                              std::array<double,CH_SPACEDIM> a_centre,
                            const int a_c_Fphi_flux, const int a_c_Sphi_source,
                             const int a_c_rho, const Interval a_c_Si,
                             const Interval a_c_Sij)
    : m_matter(a_matter), m_deriv(dx),  m_dx(dx), m_L(a_L),
                              m_centre(a_centre), m_c_Fphi_flux(a_c_Fphi_flux),
                            m_c_Sphi_source(a_c_Sphi_source), m_c_rho(a_c_rho),
                                               m_c_Si(a_c_Si), m_c_Sij(a_c_Sij)
{
    if (m_c_Si.size() != 0)
    {
        // Si is a vector
        CH_assert(m_c_Si.size() == DEFAULT_TENSOR_DIM);
    }

    if (m_c_Sij.size() != 0)
    {
        // Sij is a symmetric tensor
        CH_assert(m_c_Sij.size() ==
                  DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2);
    }
}

template <class matter_t>
template <class data_t>
void EMTensor_and_mom_flux<matter_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    if (m_c_rho >= 0)
    {
        current_cell.store_vars(emtensor.rho, m_c_rho);
    }

    if (m_c_Si.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        FOR1(i) { current_cell.store_vars(emtensor.Si[i], m_c_Si.begin() + i); }
#endif
    }

    if (m_c_Sij.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        current_cell.store_vars(emtensor.Sij[0][0], m_c_Sij.begin());
        current_cell.store_vars(emtensor.Sij[0][1], m_c_Sij.begin() + 1);
        current_cell.store_vars(emtensor.Sij[0][2], m_c_Sij.begin() + 2);
        current_cell.store_vars(emtensor.Sij[1][1], m_c_Sij.begin() + 3);
        current_cell.store_vars(emtensor.Sij[1][2], m_c_Sij.begin() + 4);
        current_cell.store_vars(emtensor.Sij[2][2], m_c_Sij.begin() + 5);
#endif
    }

    Coordinates<data_t> coords(current_cell, m_dx,m_centre);

    // start on flux term

    data_t F_phi = 0.; // angular momentum flux
    data_t aT_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // lapse times the spatial components of the 4-stress tensor
    //data_t gamma_UU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //
    data_t x = coords.x, y=coords.y, z=coords.z, r_xyz = sqrt(x*x + y*y + z*z),
                               r_xy = sqrt(x*x + y*y), sintheta = r_xy/r_xyz,
                                             sinphi = y/r_xy, cosphi = x/r_xy,
                                             costheta = z/r_xyz;
    Tensor<2, data_t, 3> gamma_UU;
    FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*vars.chi;

    FOR3(i,j,k) aT_UL[i][j] += vars.lapse*gamma_UU[i][k]*emtensor.Sij[k][j];
    FOR2(i,j) aT_UL[i][j] += -vars.shift[i]*emtensor.Si[j];

    F_phi = sintheta*(
                        cosphi*(x*aT_UL[0][1] + y*aT_UL[1][1] +z*aT_UL[2][1])
                      - sinphi*(x*aT_UL[0][0] + y*aT_UL[1][0] +z*aT_UL[2][0]))
                       /sqrt(gamma_UU[0][0]);



    // start source term

    Tensor<2, data_t, 3> J_UL;// (d x^a)/(d tildex^b)
    data_t dJ_UL[3][3][3];// tildepartial_c (d x^a)/(d tildex^b)
    data_t d_tilde_gamma_LL[3][3][3];// 3-metric 1st dervis in polars partial_c gamma_ab
    data_t tilde_gamma_LL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    data_t d_gamma_LL[3][3][3]; // cartesian gradiants of 3-metric
    data_t S_phi=0.;
    auto gamma_LL = compute_inverse_sym(gamma_UU);

    FOR3(i,j,k) dJ_UL[i][j][k] = 0.;
    FOR3(i,j,k) d_tilde_gamma_LL[i][j][k] = 0.;
    FOR3(i,j,k) d_gamma_LL[i][j][k] = d1.h[i][j][k]/vars.chi
                                      - vars.h[i][j]*d1.chi[k]*pow(vars.chi,-2);

    // checked these with mathematica
    J_UL[0][0] = cosphi*sintheta;
    J_UL[1][0] = sinphi*sintheta;
    J_UL[2][0] = costheta;
    J_UL[0][1] = r_xyz*cosphi*costheta;
    J_UL[1][1] = r_xyz*sinphi*costheta;
    J_UL[2][1] = -r_xyz*sintheta;
    J_UL[0][2] = -r_xyz*sinphi*sintheta;
    J_UL[1][2] = r_xyz*cosphi*sintheta;
    J_UL[2][2] = 0.;

    // inverse jacobean
    auto J_inv_UL = compute_inverse(J_UL); //the function outputs the transverse of inverse

    // checked these with mathematica
    // radial derivs
    dJ_UL[0][1][0] = cosphi*costheta;
    dJ_UL[1][1][0] = sinphi*costheta;
    dJ_UL[2][1][0] = -sintheta;
    dJ_UL[0][2][0] = -sinphi*sintheta;
    dJ_UL[1][2][0] = cosphi*sintheta;

    // checked these with mathematica
    // theta derivs
    dJ_UL[0][0][1] = cosphi*costheta;
    dJ_UL[1][0][1] = sinphi*costheta;
    dJ_UL[2][0][1] = -sintheta;
    dJ_UL[0][1][1] = -r_xyz*cosphi*sintheta;
    dJ_UL[1][1][1] = -r_xyz*sinphi*sintheta;
    dJ_UL[2][1][1] = -r_xyz*costheta;
    dJ_UL[0][2][1] = -r_xyz*sinphi*costheta;
    dJ_UL[1][2][1] = r_xyz*cosphi*costheta;

    // checked these with mathematica
    // phi derivs
    dJ_UL[0][0][2] = -sinphi*sintheta;
    dJ_UL[1][0][2] = cosphi*sintheta;
    dJ_UL[0][1][2] = -r_xyz*sinphi*costheta;
    dJ_UL[1][1][2] = r_xyz*cosphi*costheta;
    dJ_UL[0][2][2] = -r_xyz*cosphi*sintheta;
    dJ_UL[1][2][2] = -r_xyz*sinphi*sintheta;

    J_inv_UL[0][0] = cosphi*sintheta;
    J_inv_UL[0][1] = sinphi*sintheta;
    J_inv_UL[0][2] = costheta;
    J_inv_UL[1][0] = cosphi*costheta/r_xyz;
    J_inv_UL[1][1] = costheta*sinphi/r_xyz;
    J_inv_UL[1][2] = - sintheta/r_xyz;
    J_inv_UL[2][0] = - sinphi/(sintheta*r_xyz);
    J_inv_UL[2][1] = cosphi/(sintheta*r_xyz);
    J_inv_UL[2][2] = 0.;

    FOR4(i,j,k,l) tilde_gamma_LL[i][j] = tilde_gamma_LL[i][j] +
                                 vars.h[k][l]*(J_UL[k][i]*J_UL[l][j])/vars.chi;

    Tensor<2, data_t, 3> dummy_tensor1;
    FOR2(i,j) dummy_tensor1[i][j] = tilde_gamma_LL[i][j];

    const auto tilde_gamma_UU = compute_inverse_sym(dummy_tensor1);

    FOR5(a,b,c,n,p) d_tilde_gamma_LL[b][c][a] += 2.*dJ_UL[n][a][b]*J_UL[p][c]
                                                                *gamma_LL[n][p];

    FOR6(i,j,k,m,n,p) d_tilde_gamma_LL[i][j][k] += J_UL[m][i]*J_UL[n][j]*
                                                J_UL[p][k]*d_gamma_LL[m][n][p];

    Tensor<2, Tensor<1, data_t>> dummy_tensor2;
    FOR3(i,j,k) dummy_tensor2[i][j][k] = d_tilde_gamma_LL[i][j][k];

    const auto tilde_3_chris = compute_christoffel( dummy_tensor2, tilde_gamma_UU);

    FOR1(i) S_phi += -emtensor.rho*J_UL[i][2]*d1.lapse[i];
    FOR2(i,j) S_phi += emtensor.Si[i]*J_UL[j][2]*d1.shift[i][j];
    FOR3(i,j,m) S_phi += -J_inv_UL[i][m]*emtensor.Si[j]*dJ_UL[j][i][2]*vars.shift[m];
    FOR5(i,j,l,m,n) S_phi += vars.lapse*J_inv_UL[i][m]*J_UL[n][j]*gamma_UU[m][l]
                                  *emtensor.Sij[l][n]*tilde_3_chris.ULL[j][i][2];

    // volume elements
    data_t sqrt_tilde_gamma = sqrt(compute_determinant(dummy_tensor1)); //det of spherical polar metric
    data_t sqrt_sigma = sqrt( tilde_gamma_LL[1][1]*tilde_gamma_LL[2][2] -
                              tilde_gamma_LL[1][2]*tilde_gamma_LL[2][1] ); //det of 2 metric of surface of constant radius
    data_t sqrt_sigma_weighted = sqrt_sigma/(r_xyz*r_xy) ; // r_xyz*r_xy is r^2 sin theta


    //current_cell.store_vars(F_phi*sqrt_sigma_weighted,m_c_Fphi_flux);
    //current_cell.store_vars(S_phi*sqrt_tilde_gamma,m_c_Sphi_source); // storing S_phi * sqrt(gamma)
    current_cell.store_vars(F_phi*sqrt_sigma_weighted,m_c_Fphi_flux);
    current_cell.store_vars(S_phi*pow(vars.chi,-1.5),m_c_Sphi_source); // storing S_phi * sqrt(gamma)

}

#endif /* MOMFLUXCALC_IMPL_HPP */
