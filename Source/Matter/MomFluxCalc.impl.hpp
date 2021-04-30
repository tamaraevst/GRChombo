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
                                  const int a_c_Qphi_density, const int a_c_rho,
                                  const Interval a_c_Si, const Interval a_c_Sij)
    : m_matter(a_matter), m_deriv(dx),  m_dx(dx), m_L(a_L),
                              m_centre(a_centre), m_c_Fphi_flux(a_c_Fphi_flux),
                                              m_c_Sphi_source(a_c_Sphi_source),
                          m_c_Qphi_density(a_c_Qphi_density), m_c_rho(a_c_rho),
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

    ////////////////////////////
    // density
    ////////////////////////////

    data_t x = coords.x, y=coords.y, z=coords.z,
                               r_xyz = sqrt(x*x + y*y + z*z + 0.00000000001),
                                      r_xy = sqrt(x*x + y*y + 0.00000000001),
                                                        sintheta = r_xy/r_xyz,
                                             sinphi = y/r_xy, cosphi = x/r_xy,
                                             costheta = z/r_xyz;

    data_t Q_phi = y*emtensor.Si[0] - x*emtensor.Si[1]; // ( minus phi component of S_i)


    //////////////////////////
    // flux term
    //////////////////////////

    data_t F_phi = 0., gamma_rr_inv = 0.; // angular momentum flux
    //data_t aT_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // lapse times the spatial components of the 4-stress tensor
    data_t S_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // S^i_j

    Tensor<2, data_t, 3> gamma_UU;
    FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*vars.chi;
    auto gamma_LL = compute_inverse_sym(gamma_UU);

    //FOR3(i,j,k) aT_UL[i][j] += vars.lapse*gamma_UU[i][k]*emtensor.Sij[k][j];
    FOR3(i,j,k) S_UL[i][j] += gamma_UU[i][k]*emtensor.Sij[k][j]; // S^i_j
    //FOR2(i,j) aT_UL[i][j] += -vars.shift[i]*emtensor.Si[j];

    Tensor<2, data_t, 3> J_UL;// (d x^a)/(d tildex^b)
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

    //FOR2(i,j) gamma_rr_inv += gamma_UU[i][j]*J_inv_UL2[i][0]*J_inv_UL2[j][0]; same as explicit version, checked

    gamma_rr_inv = (x*x*gamma_UU[0][0] + y*y*gamma_UU[1][1] + z*z*gamma_UU[2][2]
               + 2.*(x*y*gamma_UU[0][1] + y*z*gamma_UU[1][2] + x*z*gamma_UU[0][2])
             )/(r_xyz*r_xyz); // same as version above, checked. This is gamma^{rr}
             
    /*F_phi = sintheta*(
                        cosphi*(x*aT_UL[0][1] + y*aT_UL[1][1] +z*aT_UL[2][1])
                      - sinphi*(x*aT_UL[0][0] + y*aT_UL[1][0] +z*aT_UL[2][0])
                    )/sqrt(gamma_rr_inv);*/ //just a different way of expressing it

    /*F_phi = (   x*(x*aT_UL[0][1] + y*aT_UL[1][1] +z*aT_UL[2][1])
              - y*(x*aT_UL[0][0] + y*aT_UL[1][0] +z*aT_UL[2][0])
            )/(sqrt(gamma_rr_inv)*r_xyz);*/


    F_phi = (  x*( x*S_UL[0][1] + y*S_UL[1][1] + z*S_UL[2][1] )
              - y*( x*S_UL[0][0] + y*S_UL[1][0] + z*S_UL[2][0] )
                                        )*vars.lapse/(sqrt(gamma_rr_inv)*r_xyz);
    F_phi -= ( x*vars.shift[0] + y*vars.shift[1] + z*vars.shift[2] )
                *(x*emtensor.Si[1]-y*emtensor.Si[0])/(sqrt(gamma_rr_inv)*r_xyz);

    /*data_t kroneka[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    FOR3(i,j,k) kroneka[i][k] += J_UL2[i][j]*J_inv_UL2[k][j];*/ // checked inverse metric works



    /////////////////////////////////////
    // start source term again in cartesian
    /////////////////////////////////////

    data_t S_phi = 0.;
    Tensor<1, data_t, 3> xi; // approx killing vector = partial_phi but expressed in cartesian coords
    xi[0] = -y;
    xi[1] = x;
    xi[2]=0.;
    Tensor<2, data_t, 3> d_xi;
    FOR2(i,j) d_xi[i][j]=0.; //deriv is i, xi is j (opposite to other deriv convention for no reason) // partial_i xi^j
    d_xi[1][0] = -1.;
    d_xi[0][1] = 1.;

    FOR1(i) S_phi += - emtensor.rho*xi[i]*d1.lapse[i];
    FOR2(i,j) S_phi += emtensor.Si[i]*(xi[j]*d1.shift[i][j] - vars.shift[j]*d_xi[j][i]);  // assumed d1.shift[i][j] = partial_j beta^i
    FOR3(i,j,k) S_phi += vars.lapse*emtensor.Sij[i][j]*gamma_UU[j][k]*d_xi[k][i];
    FOR4(i,j,k,l) S_phi += vars.lapse*emtensor.Sij[j][k]*gamma_UU[k][i]*chris.ULL[j][i][l]*xi[l];



    /////////////////////////////////////
    // volumes/areas
    /////////////////////////////////////
    data_t polar_metric[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    FOR4(i,j,k,l) polar_metric[i][j] += gamma_LL[k][l]*J_UL[k][i]*J_UL[l][j]; // make polar 3-metric

    data_t sqrt_sigma = sqrt( polar_metric[1][1]*polar_metric[2][2] -
                              polar_metric[1][2]*polar_metric[2][1] ); //det of 2 metric of surface of constant radius

    data_t sqrt_gamma_polar = sqrt(
      polar_metric[0][0]*( polar_metric[1][1]*polar_metric[2][2] - polar_metric[1][2]*polar_metric[2][1] ) +
      polar_metric[0][1]*( polar_metric[1][2]*polar_metric[2][0] - polar_metric[2][2]*polar_metric[1][0] ) +
      polar_metric[0][2]*( polar_metric[1][0]*polar_metric[2][1] - polar_metric[1][1]*polar_metric[2][0] )
    ); //det of polar 3-metric

    data_t alternative_sqrt_sigma = sqrt(gamma_rr_inv*sqrt_gamma_polar*sqrt_gamma_polar); // dont need this, but checked its equal to sqrt_sigma

    data_t sqrt_gamma = pow(vars.chi,-1.5); //det of cartesian 3-metric

    data_t sqrt_sigma_modified = sqrt_sigma/(r_xyz*r_xy); // r_xyz*r_xy is r^2 sin theta that is automatically in flux int

    current_cell.store_vars(Q_phi*sqrt_gamma,m_c_Qphi_density);
    current_cell.store_vars(F_phi*sqrt_sigma_modified,m_c_Fphi_flux); // rememebr, the spherical integrator that uses this already includes r^2*sin theta
    current_cell.store_vars(S_phi*sqrt_gamma,m_c_Sphi_source); // storing S_phi * sqrt(gamma)


}

#endif /* MOMFLUXCALC_IMPL_HPP */
