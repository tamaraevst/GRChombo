//Last update K Clough 23.04.2017

#if !defined(KERRBH_HPP_)
#error "This file should only be included through KerrBH.hpp"
#endif

#ifndef KERRBH_IMPL_HPP_
#define KERRBH_IMPL_HPP_

// Computes semi-isotropic Kerr solution as detailed in Liu, Etienne and Shapiro 2010, 
// arxiv gr-qc/1001.4077
template <class data_t>
void KerrBH::compute(Cell current_cell)
{
    // set up vars for the metric and extrinsic curvature, shift and lapse in spherical coords
    tensor<2,data_t> spherical_g;
    tensor<2,data_t> spherical_K;
    tensor<1,data_t> spherical_shift;
    data_t kerr_lapse;

    // The cartesian variables and coords
    Vars<data_t> vars;
    Coordinates<data_t> coords(current_cell, m_dx);

    // Compute the components in spherical coords as per 1401.1548
    compute_kerr(spherical_g, spherical_K, spherical_shift, kerr_lapse, coords);

    // get position
    data_t x, r;
    double y, z;
    get_position(coords, x, y, z, r);

    using namespace InitialDataTools;
    // Convert spherical components to cartesian components using coordinate transforms
    vars.h = spherical_to_cartesian_LL(spherical_g, r, x, y, z);
    vars.A = spherical_to_cartesian_LL(spherical_K, r, x, y, z);
    vars.shift = spherical_to_cartesian_U(spherical_shift, r, x, y, z);

    using namespace TensorAlgebra;
    //Convert to BSSN vars
    data_t deth = compute_determinant(vars.h);
    auto h_UU = compute_inverse(vars.h);
    vars.chi = pow(deth, -1./3.);

    // transform extrinsic curvature into A and TrK - note h is still non conformal version
    // which is what we need here
    vars.K = compute_trace(vars.A, h_UU);
    make_trace_free(vars.A, vars.h, h_UU);

    //Make conformal
    FOR2(i,j)
    {
        vars.h[i][j] *= vars.chi;
        vars.A[i][j] *= vars.chi;
    }

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
    get_position(coords, x, y, z, r);
    r2 = r*r;

    //the radius in xy plane, subject to a floor
    rho2 = pow(x, 2.0) + pow(y, 2.0);
    double minimum_rho2 = 1e-12;
    auto rho_is_too_small = simd_compare_lt(rho2, minimum_rho2);
    rho2 = simd_conditional(rho_is_too_small, minimum_rho2, rho2);
    rho = sqrt(rho2);

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
void KerrBH::get_position(Coordinates<data_t> coords, data_t &x, double &y, double &z, data_t &r)
{
    // work out where we are on the grid
    x = coords.x - m_params.center[0];
    y = coords.y - m_params.center[1];
    z = coords.z - m_params.center[2];

    //the radius, subject to a floor
    r = coords.get_radius(m_params.center);
}

#endif /* KERRBH_IMPL_HPP_ */
