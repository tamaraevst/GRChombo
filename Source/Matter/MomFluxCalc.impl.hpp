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

    const auto phys_chris = compute_phys_chris(d1.chi, vars.chi, vars.h,
                                                          h_UU, chris.ULL);

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
    // coordinates and killing vector
    ////////////////////////////

    data_t x = coords.x, y=coords.y, z=coords.z,
                               r_xyz = sqrt(x*x + y*y + z*z + 0.00000000001),
                                      r_xy = sqrt(x*x + y*y + 0.00000000001),
                                                        sintheta = r_xy/r_xyz,
                                             sinphi = y/r_xy, cosphi = x/r_xy,
                                             costheta = z/r_xyz;
    Tensor<1, data_t, 3> xi; // approx killing vector = partial_phi but expressed in cartesian coords
    xi[0] = -y;
    xi[1] = x;
    xi[2]=0.;
    Tensor<1, data_t, 3> cart_coords;
    cart_coords[0] = x;
    cart_coords[1] = y;
    cart_coords[2] = z;
    Tensor<2, data_t, 3> d_xi;
    FOR2(i,j) d_xi[i][j]=0.; //deriv is i, xi is j (opposite to other deriv convention for no reason) // partial_i xi^j
    d_xi[1][0] = -1.;
    d_xi[0][1] = 1.;







    ////////////////////////////
    // metrics, volumes
    ////////////////////////////

    Tensor<2, data_t, 3> gamma_UU;
    FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*vars.chi;
    auto gamma_LL = compute_inverse_sym(gamma_UU);

    Tensor<2, data_t, 3> J_UL;// (d x^a)/(d tildex^b) jacobean
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

    /*data_t kroneka[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    FOR3(i,j,k) kroneka[i][k] += J_UL2[i][j]*J_inv_UL2[k][j];*/ // checked inverse metric works

    data_t gamma_polar_LL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // downstairs indeces
    Tensor<2, data_t, 3> dummy; //placeholder equal to polar metric
    FOR4(i,j,k,l) gamma_polar_LL[i][j] += gamma_LL[k][l]*J_UL[k][i]*J_UL[l][j]; // make spatial polar 3-metric
    FOR2(i,j) dummy[i][j] = gamma_polar_LL[i][j];
    auto gamma_polar_UU = compute_inverse_sym(dummy);

    data_t sqrt_sigma = sqrt( gamma_polar_LL[1][1]*gamma_polar_LL[2][2] -
                              gamma_polar_LL[1][2]*gamma_polar_LL[2][1] ); //det of 2 metric of surface of constant radius and time

    data_t sqrt_sigma_modified = sqrt_sigma/(r_xyz*r_xy); // r_xyz*r_xy is r^2 sin theta that is automatically in flux int

    // cylinder metric
    Tensor<2, data_t, 3> hamma_LL; // cant use h, so called it hamma. its the 3-cylinder metric
    Tensor<1, data_t, 3> beta_L; //cartesian
    Tensor<1, data_t, 3> beta_L_polar;
    FOR1(i) beta_L[i]=0.;
    FOR1(i) beta_L_polar[i]=0.;
    FOR2(i,j) beta_L[i] += gamma_LL[i][j]*vars.shift[j];
    FOR2(i,j) beta_L_polar[i] += J_UL[j][i]*beta_L[j];
    hamma_LL[0][0] = -vars.lapse*vars.lapse;
    FOR2(i,j) hamma_LL[0][0] += vars.shift[i]*vars.shift[j]*gamma_LL[i][j];
    hamma_LL[1][1] = gamma_polar_LL[1][1];
    hamma_LL[1][2] = gamma_polar_LL[1][2];
    hamma_LL[2][1] = gamma_polar_LL[2][1];
    hamma_LL[2][2] = gamma_polar_LL[2][2];
    hamma_LL[1][0] = beta_L_polar[1];
    hamma_LL[2][0] = beta_L_polar[2];
    hamma_LL[0][1] = beta_L_polar[1];
    hamma_LL[0][2] = beta_L_polar[2];
    auto hamma_UU = compute_inverse_sym(hamma_LL);
    data_t root_minus_h = sqrt(-compute_determinant_sym(hamma_LL));
    data_t root_minus_h_modified = root_minus_h/(r_xyz*r_xy); // modified by 1/(r^2 sintheta) which is in the flux integrator already
    //data_t hopefully_zero = 1. - sqrt_sigma/(root_minus_h*sqrt(-hamma_UU[0][0])); // it is :)

    data_t beta_r = ( x*vars.shift[0]+ y*vars.shift[1]+ z*vars.shift[2]  )/r_xyz; //upstars component
    data_t g_rr_inv = gamma_polar_UU[0][0] - beta_r*beta_r/(vars.lapse*vars.lapse);
    //data_t hopefully_zero = vars.lapse/sqrt(gamma_rr_inv) - 1./sqrt(-hamma_UU[0][0]*g_rr_inv); // it is :)
    //data_t hopefully_zero = 1. + gamma_rr_inv/(hamma_UU[0][0]*g_rr_inv*vars.lapse*vars.lapse); // it is :)


    data_t sqrt_gamma_polar = sqrt(
      gamma_polar_LL[0][0]*( gamma_polar_LL[1][1]*gamma_polar_LL[2][2] - gamma_polar_LL[1][2]*gamma_polar_LL[2][1] ) +
      gamma_polar_LL[0][1]*( gamma_polar_LL[1][2]*gamma_polar_LL[2][0] - gamma_polar_LL[2][2]*gamma_polar_LL[1][0] ) +
      gamma_polar_LL[0][2]*( gamma_polar_LL[1][0]*gamma_polar_LL[2][1] - gamma_polar_LL[1][1]*gamma_polar_LL[2][0] )
    ); //det of polar 3-metric
    //data_t alternative_sqrt_sigma = sqrt(gamma_rr_inv*sqrt_gamma_polar*sqrt_gamma_polar); // dont need this, but checked its equal to sqrt_sigma

    //data_t sqrt(compute_determinant_sym(gamma_LL)); // same as line below
    data_t sqrt_gamma = pow(vars.chi,-1.5); //det of cartesian 3-metric








    ////////////////////////////
    // density
    ////////////////////////////

    data_t Q_phi = 0.; //phi component of em tensors Si
    FOR1(i) Q_phi += -xi[i]*emtensor.Si[i];
    //data_t Q_phi = y*emtensor.Si[0] - x*emtensor.Si[1]; // ( minus phi component of S_i)






    //////////////////////////
    // flux term
    //////////////////////////

    // the polar version with sqrt(sigma) as det
    // this version needs int F_phi sqrt(sigma) dx^2
    // dont forget to divide sqrt(sigma) by r_xy*r_xyz = r^2 sintheta as the surface integrator has r^2 sin theta already

    data_t F_phi = 0.;
    data_t S_UL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // S^i_j
    FOR3(i,j,k) S_UL[i][j] += gamma_UU[i][k]*emtensor.Sij[k][j]; // S^i_j
    data_t S_mixed[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; // S_ij with i polar j cart
    FOR3(i,j,k) S_mixed[i][j] += emtensor.Sij[k][j]*J_UL[k][j];

    //FOR2(i,j) F_phi += vars.lapse*(cart_coords[i]/r_xyz)*S_UL[i][j]*xi[j];
    //FOR2(i,j) F_phi += vars.lapse*gamma_polar_UU[0][i]*S_mixed[i][j]*xi[j];
    FOR1(j) F_phi += vars.lapse*S_mixed[0][j]*xi[j]/sqrt(gamma_polar_LL[0][0]); //why is this LL, shoudl be UU?
    FOR1(i) F_phi += - beta_r*emtensor.Si[i]*xi[i]/sqrt(gamma_polar_UU[0][0]);
    //F_phi = F_phi/sqrt(gamma_polar_UU[0][0]);









    //////////////////////////
    // flux term again
    //////////////////////////

    // the cartesian version with sqrt(-h) as det
    // this version requires int F sqrt(-h) dx^2
    // dont forget to divide sqrt(-h) by r_xy*r_xyz = r^2 sintheta as this term is automatically in the surface integrator

    data_t F2_phi = 0., N_normsqr = 0.; // angular momentum flux
    Tensor<2, data_t, 3> g_UU; //4-metric but spatial parts only (not this is not a proper metric)
    FOR2(i,j) g_UU  = gamma_UU[i][j] - vars.shift[i]*vars.shift[j]/(vars.lapse*vars.lapse);
    FOR2(i,j) N_normsqr += g_UU[i][j]*cart_coords[i]*cart_coords[j];
    Tensor<1, data_t, 3> N_cart; // downstaits componnet of radial unit vector
    Tensor<1, data_t, 3> N_proj; // upstairs componetnets projected to Sigma
    N_cart[0] = x/sqrt(N_normsqr);
    N_cart[1] = y/sqrt(N_normsqr);
    N_cart[2] = z/sqrt(N_normsqr);
    N_proj[0] = 0.;
    N_proj[1] = 0.;
    N_proj[2] = 0.;
    FOR2(i,j) N_proj[i] += gamma_UU[i][j]*N_cart[j];

    FOR2(i,j) F2_phi += -N_cart[i]*vars.shift[i]*xi[j]*emtensor.Si[j]/vars.lapse;
    FOR2(i,j) F2_phi += N_proj[i]*xi[j]*emtensor.Sij[i][j];




    //////////////////////////
    // flux term againnnnn, this one seems to work!
    //////////////////////////

    data_t F3_phi=0.;





    data_t Sij_polar[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    data_t Si_polar[3] = {0.,0.,0.};
    FOR2(i,j) Si_polar[i] += emtensor.Si[j]*J_UL[j][i];
    FOR4(i,j,k,l) Sij_polar[i][j] += emtensor.Sij[k][l]*J_UL[k][i]*J_UL[l][j];
    FOR1(i) F3_phi += vars.lapse*Sij_polar[i][2]*gamma_polar_UU[0][i];
    F3_phi += - beta_r*Si_polar[2];
    F3_phi = F3_phi/sqrt(gamma_polar_UU[0][0]);




    /////////////////////////////////////
    // start source in cartesian
    /////////////////////////////////////

    data_t S_phi = 0.;

    FOR1(i) S_phi += - emtensor.rho*xi[i]*d1.lapse[i];
    FOR2(i,j) S_phi += emtensor.Si[i]*(xi[j]*d1.shift[i][j] - vars.shift[j]*d_xi[j][i]);  // assumed d1.shift[i][j] = partial_j beta^i
    FOR3(i,j,k) S_phi += vars.lapse*emtensor.Sij[i][j]*gamma_UU[j][k]*d_xi[k][i];
    FOR4(i,j,k,l) S_phi += vars.lapse*emtensor.Sij[j][k]*gamma_UU[k][i]*chris.ULL[j][i][l]*xi[l];



    /////////////////////////////////////
    // Source term again
    /////////////////////////////////////
    data_t S2_phi=0.;
    data_t XI[3] = {-y,x,0.};
    data_t DXI[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    data_t CHRIS3[3][3][3];
    DXI[0][1] = 1.;
    DXI[1][0] = -1.;
    FOR1(i) S2_phi += -emtensor.rho*XI[i]*d1.lapse[i];
    FOR2(i,j) S2_phi += emtensor.Si[i]*XI[j]*d1.shift[i][j] - emtensor.Si[j]*vars.shift[i]*DXI[i][j];
    FOR3(i,j,k) S2_phi += vars.lapse*gamma_UU[i][j]*emtensor.Sij[j][k]*DXI[i][k];
    FOR4(i,j,k,l) S2_phi += vars.lapse*gamma_UU[i][j]*emtensor.Sij[j][k]*phys_chris[k][i][l]*XI[l];







    /////////////////////////////////////
    // store values
    /////////////////////////////////////

    current_cell.store_vars(Q_phi*sqrt_gamma,m_c_Qphi_density);
    //current_cell.store_vars(F2_phi*root_minus_h_modified - F_phi*sqrt_sigma_modified,m_c_Fphi_flux); // rememebr, the spherical integrator that uses this already includes r^2*sin theta
    current_cell.store_vars(F3_phi*sqrt_sigma_modified,m_c_Fphi_flux); // rememebr, the spherical integrator that uses this already includes r^2*sin theta
    //current_cell.store_vars(is_this_zero,m_c_Fphi_flux); // rememebr, the spherical integrator that uses this already includes r^2*sin theta
    current_cell.store_vars(S2_phi*sqrt_gamma,m_c_Sphi_source); // storing S_phi * sqrt(gamma)






}

#endif /* MOMFLUXCALC_IMPL_HPP */
