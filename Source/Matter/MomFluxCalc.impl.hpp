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
                              std::array<double,CH_SPACEDIM> a_centre, const int a_c_Fx_flux,
                             const int a_c_Fy_flux, const int a_c_Sx_source,
                             const int a_c_Sy_source, const int a_c_rho,
                             const Interval a_c_Si, const Interval a_c_Sij)
    : m_matter(a_matter), m_deriv(dx),  m_dx(dx), m_L(a_L),
                                  m_centre(a_centre), m_c_Fx_flux(a_c_Fx_flux),
                        m_c_Fy_flux(a_c_Fy_flux), m_c_Sx_source(a_c_Sx_source),
                                m_c_Sy_source(a_c_Sy_source), m_c_rho(a_c_rho),
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

    //double N[3] = {coords.x,coords.y,coords.z};
    /*data_t N[3] = {coords.x,coords.y,coords.z};
    data_t normsqr = 0.;
    data_t Fx=0., Fy=0., Sx =0., Sy = 0.;

    FOR2(i,j) normsqr += N[i]*N[j]*vars.h[i][j]/vars.chi;
    FOR1(i) N[i] /= sqrt(normsqr);
    FOR1(i) Fx += vars.lapse*N[i]*emtensor.Sij[i][0];
    FOR1(i) Fy += vars.lapse*N[i]*emtensor.Sij[i][1];
    FOR2(i,j) Fx -=  vars.h[i][j]*vars.shift[i]*N[j]*emtensor.Si[0]/vars.chi;
    FOR2(i,j) Fy -=  vars.h[i][j]*vars.shift[i]*N[j]*emtensor.Si[1]/vars.chi;*/

    // start on flux term

    data_t azimuthal_radial_shear = 0.; // the T^r_phi component
    data_t aT_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // lapse times the spatial components of the 4-stress tensor
    data_t gamma_UU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //
    data_t x = coords.x, y=coords.y, z=coords.z, r_xyz = sqrt(x*x + y*y + z*z),
                               r_xy = sqrt(x*x + y*y), sintheta = r_xy/r_xyz,
                                             sinphi = y/r_xy, cosphi = x/r_xy,
                                             costheta = z/r_xyz;

    FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*vars.chi;

    FOR3(i,j,k) aT_UL[i][j] += vars.lapse*gamma_UU[i][k]*emtensor.Sij[k][j];
    FOR2(i,j) aT_UL[i][j] += -vars.shift[i]*emtensor.Si[j];

    azimuthal_radial_shear = sintheta*(
                        cosphi*(x*aT_UL[0][1] + y*aT_UL[1][1] +z*aT_UL[2][1])
                      - sinphi*(x*aT_UL[0][0] + y*aT_UL[1][0] +z*aT_UL[2][0]));



    // start source term

    data_t J_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};// (d x^a)/(d tildex^b)
    data_t dJ_UL[3][3][3];// tildepartial_c (d x^a)/(d tildex^b)
    data_t d_tilde_gamma_LL[3][3][3];// 3-metric 1st dervis in polars partial_c gamma_ab
    data_t tilde_gamma_LL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    data_t S_phi=0.;

    FOR3(i,j,k) dJ_UL[i][j][k] = 0.;
    FOR3(i,j,k) d_tilde_gamma_LL[i][j][k] = 0.;

    J_UL[0][0] = cosphi*sintheta;
    J_UL[1][0] = sinphi*sintheta;
    J_UL[2][0] = costheta;
    J_UL[0][1] = r_xyz*cosphi*costheta;
    J_UL[1][1] = r_xyz*sinphi*costheta;
    J_UL[2][1] = -r_xyz*sintheta;
    J_UL[0][2] = -r_xyz*sinphi*sintheta;
    J_UL[1][2] = r_xyz*cosphi*sintheta;

    // radial derivs
    dJ_UL[0][1][0] = cosphi*costheta;
    dJ_UL[1][1][0] = sinphi*costheta;
    dJ_UL[2][1][0] = -sintheta;
    dJ_UL[0][2][0] = -sinphi*sintheta;
    dJ_UL[1][2][0] = cosphi*sintheta;

    //theta derivs
    dJ_UL[0][0][1] = cosphi*costheta;
    dJ_UL[1][0][1] = sinphi*costheta;
    dJ_UL[2][0][1] = -sintheta;
    dJ_UL[0][1][1] = -r_xyz*cosphi*sintheta;
    dJ_UL[1][1][1] = -r_xyz*sinphi*sintheta;
    dJ_UL[2][1][1] = -r_xyz*costheta;
    dJ_UL[0][2][1] = -r_xyz*sinphi*costheta;
    dJ_UL[1][2][1] = r_xyz*cosphi*costheta;

    //phi derivs
    dJ_UL[0][0][2] = -sinphi*sintheta;
    dJ_UL[1][0][2] = cosphi*sintheta;
    dJ_UL[0][1][2] = -r_xyz*sinphi*costheta;
    dJ_UL[1][1][2] = r_xyz*cosphi*costheta;
    dJ_UL[0][2][2] = -r_xyz*cosphi*sintheta;
    dJ_UL[1][2][2] = -r_xyz*sinphi*sintheta;

    FOR4(i,j,k,l) tilde_gamma_LL[i][j] += vars.h[k][l]*
                                               (J_UL[k][i]*J_UL[l][j])/vars.chi;

    auto tilde_gamma_UU = compute_inverse_sym(tilde_gamma_LL);




    current_cell.store_vars(azimuthal_radial_shear,m_c_Fx_flux);
    current_cell.store_vars(0.,m_c_Fy_flux);
    current_cell.store_vars(S_phi,m_c_Sx_source);
    current_cell.store_vars(0.,m_c_Sy_source);

}

#endif /* MOMFLUXCALC_IMPL_HPP */
