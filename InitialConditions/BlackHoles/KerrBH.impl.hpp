//Last update K Clough 23.04.2017

#if !defined(KERRBH_HPP_)
#error "This file should only be included through KerrBH.hpp"
#endif

#ifndef KERRBH_IMPL_HPP_
#define KERRBH_IMPL_HPP_

#include "KerrBH.hpp"
#include "simd.hpp"
#include "DebuggingTools.hpp"

// Computes semi-isotropic Kerr solution as detailed in Liu, Etienne and Shapiro 2010, 
// arxiv gr-qc/1001.4077
template <class data_t>
void KerrBH::compute(Cell current_cell)
{
    // set up vars for the metric and extrinsic curvature in spherical coords
    tensor<2,data_t> spherical_g;
    tensor<2,data_t> spherical_K;

    // this is the spherical shift (phi is the only non zero one in Kerr, but function is general)
    tensor<1,data_t> spherical_shift;

    // the analytic lapse
    data_t kerr_lapse;

    // The other variables
    Vars<data_t> vars;

    // Create a Coordinates object for current cell
    Coordinates<data_t> coords(current_cell, m_dx);

    // Compute the components in spherical coords as per 1401.1548
    compute_kerr(spherical_g, spherical_K, spherical_shift, kerr_lapse, coords);

    // Convert spherical components to cartesian components using coordinate transforms
    spherical_to_cartesian(vars, spherical_g, spherical_K, spherical_shift, coords);

    // use a pre collapsed lapse, could also use analytic one
    // vars.lapse = kerr_lapse;
    vars.lapse = pow(vars.chi, 0.5);

    // Populate the variables on the grid
    // NB We stil need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function 
    // as we need the gradients of the metric which are not yet available
    current_cell.store_vars(vars);
}

template <class data_t>
void KerrBH::compute_kerr(tensor<2,data_t> &spherical_g,
                          tensor<2,data_t> &spherical_K,
                          tensor<1,data_t> &spherical_shift,
                          data_t &kerr_lapse,
                          const Coordinates<data_t> coords)
{   
    //Kerr black hole params - mass M and spin a
    double M = m_params.mass;
    double a = m_params.spin;

    // get position
    data_t x, r, r2, rho, rho2;
    double y, z;
    get_position(coords, x, y, z, r, r2, rho, rho2);

    //calculate useful position quantities
    data_t costh = z/r;
    data_t sinth = rho/r;
    data_t costh2 = costh*costh;
    data_t sinth2 = sinth*sinth;

    // calculate useful metric quantities
    double r_plus = M + sqrt(M*M - a*a);
    double r_minus = M - sqrt(M*M - a*a);

    //The Boyer-Lindquist coordinate
    data_t r_BL = r*pow(1.0 + 0.25*r_plus/r, 2.0);

    //Other useful quantities per 1001.4077
    data_t Sigma = r_BL*r_BL + a*a*costh2;
    data_t Delta = r_BL*r_BL - 2.0*M*r_BL + a*a;
    // In the paper this is just 'A', but not to be confused with A_ij
    data_t AA = pow(r_BL*r_BL + a*a, 2.0) - Delta*a*a*sinth2;
    // The rr component of the conformal spatial matric
    data_t gamma_rr = Sigma*pow(r + 0.25*r_plus, 2.0)/(r*r2*(r_BL - r_minus));

    // Metric in semi isotropic Kerr-Schild coordinates, r, theta (t or th), phi (p)
    FOR2(i,j)
    {
        spherical_g[i][j] = 0.0;
    }
    spherical_g[0][0] = gamma_rr;           // gamma_rr
    spherical_g[1][1] = Sigma;              // gamma_tt
    spherical_g[2][2] = AA/Sigma * sinth2;  // gamma_pp

    // Extrinsic curvature
    FOR2(i,j)
    {
        spherical_K[i][j] = 0.0;
    }

    // set non zero elements of Krtp - K_rp, K_tp
    spherical_K[0][2] = a*M*sinth2/(Sigma*sqrt(AA*Sigma))*
                          (3.0*pow(r_BL, 4.0) + 2*a*a*r_BL*r_BL - pow(a,4.0) - a*a*(r_BL*r_BL - a*a)*sinth2)
                                                         *(1.0+0.25*r_plus/r)/sqrt(r*r_BL - r*r_minus);
    spherical_K[2][0] = spherical_K[0][2];
    spherical_K[2][1] = -2.0*pow(a,3.0)*M*r_BL*costh*sinth*sinth2/(Sigma*sqrt(AA*Sigma)) 
                                                         *(r - 0.25*r_plus)*sqrt(r_BL/r - r_minus/r);
    spherical_K[1][2] = spherical_K[2][1];

    // set the analytic lapse
    kerr_lapse = sqrt(Delta*Sigma/AA);

    // set the shift (only the phi component is non zero)
    spherical_shift[0] = 0.0;
    spherical_shift[1] = 0.0;
    spherical_shift[2] = - 2.0*M*a*r_BL / AA;
}

template <class data_t>
void KerrBH::spherical_to_cartesian(Vars<data_t> &vars,
                                    tensor<2,data_t> &spherical_g,
                                    tensor<2,data_t> &spherical_K,
                                    tensor<1,data_t> &spherical_shift,
                                    const Coordinates<data_t> coords)
{
    // get position
    data_t x, r, r2, rho, rho2;
    double y, z;
    get_position(coords, x, y, z, r, r2, rho, rho2);

    // calculate useful position quantities
    data_t costh = z/r;
    data_t sinth = rho/r;
    data_t costh2 = costh*costh;
    data_t sinth2 = sinth*sinth;
    data_t cosph = x/rho;
    data_t sinph = y/rho;
    data_t cosph2 = cosph*cosph;
    data_t sinph2 = sinph*sinph;

    // derivatives for Jacobian matrix - drdx etc
    tensor<2,data_t> jac;
    jac[0][0] = x/r;
    jac[1][0] = cosph*z/r2;
    jac[2][0] = -y/rho2;
    jac[0][1] = y/r;
    jac[1][1] = sinph*z/r2;
    jac[2][1] = x/rho2;
    jac[0][2] = z/r;
    jac[1][2] = -rho/r2;
    jac[2][2] = 0.0;

    // work out detj for the jacobian, needed to transform tensor densities
    data_t detj = jac[0][0]*(jac[1][1]*jac[2][2]-jac[1][2]*jac[2][1])-
                  jac[0][1]*(jac[2][2]*jac[1][0]-jac[1][2]*jac[2][0])+
                  jac[0][2]*(jac[1][0]*jac[2][1]-jac[1][1]*jac[2][0]);

    // make conformal before transforming (could also do after)
    data_t detg = spherical_g[0][0]*(spherical_g[1][1]*spherical_g[2][2]-spherical_g[1][2]*spherical_g[2][1])-
                  spherical_g[0][1]*(spherical_g[2][2]*spherical_g[1][0]-spherical_g[1][2]*spherical_g[2][0])+
                  spherical_g[0][2]*(spherical_g[1][0]*spherical_g[2][1]-spherical_g[1][1]*spherical_g[2][0]);

    FOR2(i,j)
    { 
        spherical_g[i][j] = spherical_g[i][j]*pow(detg, -1./3.);
    }

    // Set all components to zero
    vars.assign(0.);

    // now transform using detj to allow for tensor densities
    vars.chi = pow(detg*detj*detj, -1./3.);

    // transform metric into xyz coordinates
    FOR4(i,j,k,m)
    {
        vars.h[i][j] += spherical_g[k][m]*jac[k][i]*jac[m][j]*pow(detj, -2./3.);
    }

    // invert conformal metric
    using namespace TensorAlgebra;
    auto h_UU = compute_inverse(vars.h);

    // transform extrinsic curvature into xyz coordinates and calculate A and TrK
    FOR4(i,j,k,m)
    {
        vars.A[i][j] += spherical_K[k][m]*jac[k][i]*jac[m][j];
    }

    // now make into transverse traceless form
    FOR2(i,j)
    {
        vars.K += vars.A[i][j]*h_UU[i][j]*vars.chi;
    }

    FOR2(i,j)
    {
        vars.A[i][j] = vars.A[i][j]*vars.chi - vars.h[i][j]*vars.K/3.0;
    }

    // The shift - NB a vector not a one form so need inverse - dxdphi, dydphi etc instead
    tensor<2,data_t> inv_jac;
    inv_jac[0][0] = x/r;
    inv_jac[1][0] = y/r;
    inv_jac[2][0] = z/r;
    inv_jac[0][1] = z*cosph;
    inv_jac[1][1] = z*sinph;
    inv_jac[2][1] = -r*sinth;
    inv_jac[0][2] = -y;
    inv_jac[1][2] = x;
    inv_jac[2][2] = 0.0;

    // transform the shift
    FOR1(i)
    {
        vars.shift[i] = 0.0;
        FOR1(j)
        {
            vars.shift[i] += inv_jac[i][j]*spherical_shift[j];
        }
    }
}

template <class data_t>
void KerrBH::get_position(Coordinates<data_t> coords, data_t &x, double &y, double &z, 
                          data_t &r, data_t &r2, data_t &rho, data_t &rho2)
{
    // work out where we are on the grid
    x = coords.x - m_params.center[0];
    y = coords.y - m_params.center[1];
    z = coords.z - m_params.center[2];

    //the radius, subject to a floor
    r = coords.get_radius(m_params.center);
    r2 = r*r;

    //the radius in xy plane, subject to a floor
    rho2 = pow(x, 2.0) + pow(y, 2.0);
    double minimum_rho2 = 1e-12;
    auto rho_is_too_small = simd_compare_lt(rho2, minimum_rho2);
    rho2 = simd_conditional(rho_is_too_small, minimum_rho2, rho2);
    rho = sqrt(rho2);
}

#endif /* KERRBH_IMPL_HPP_ */
