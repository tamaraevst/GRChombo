/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ROTATINGBOSONSTAR_HPP_)
#error "This file should only be included through RotatingBosonStar.hpp"
#endif

#ifndef ROTATINGBOSONSTAR_IMPL_HPP_
#define ROTATINGBOSONSTAR_IMPL_HPP_

#include <cmath>
#include "RotatingBosonStarSolution.hpp"
#include "DebuggingTools.hpp"

inline RotatingBosonStar::RotatingBosonStar(RotatingBosonStar_params_t a_params_RotatingBosonStar, double a_dx)
    : m_params_RotatingBosonStar(a_params_RotatingBosonStar), m_dx(a_dx)
{
}

void RotatingBosonStar::compute_1d_rotating_solution()
{   
    try
    {  
        rotating_BS_sol.main(m_params_RotatingBosonStar.base_path);
        pout() << "Wooo I have read the spinning BS initial data!" << endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void RotatingBosonStar::compute(Cell<data_t> current_cell) const
{   
    double theta, phi;

    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_RotatingBosonStar.star_centre);

    // Star positioning
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    theta = acos(z/r);
    phi = atan2(y,x);

    // Compactified coordinate
    double xvar = r / (1. + r);

    double ff = m_params_RotatingBosonStar.BS_frequency;
    // rotating_BS_sol.get_BSfrequency();

    int n = rotating_BS_sol.n;
    int m = rotating_BS_sol.m;
    
    double A_val = rotating_BS_sol.get_amp_interp(xvar, theta, n, m)*sqrt(2);
    // DEBUG_OUT(A_val);
    double f_val = rotating_BS_sol.get_f_interp(xvar, theta, n, m);
    double g_val = rotating_BS_sol.get_g_interp(xvar, theta, n, m);
    double l_val = rotating_BS_sol.get_l_interp(xvar, theta, n, m);
    double omega_val = rotating_BS_sol.get_omega_interp(xvar, theta, n, m);
    double dthomega_val = rotating_BS_sol.get_dthomega_interp(xvar, theta, n, m);
    double dromega_val = rotating_BS_sol.get_dromega_interp(xvar, theta, n, m);

    double lapse = sqrt(fabs(f_val));
    double beta_x = sin(theta)*sin(phi)*omega_val;
    double beta_y = -sin(theta)*cos(phi)*omega_val;
    double beta_z = 0;
    // double beta_phi = -(l_val/f_val) * r * omega_val * sin(theta) * sin(theta); 
    double beta_phi = - omega_val/r; 

    vars.shift[0] += beta_x;
    vars.shift[1] += beta_y;
    vars.shift[2] +- beta_z;

    // pout() << "Computed lapse and shift expressions \n" << endl;

    double phase_ = phi;

    double g_zz_1 = (g_val * l_val) / f_val;
    double g_yy_1 = l_val / f_val * (cos(phi)*cos(phi) + sin(phi)*sin(phi)*g_val);
    double g_xx_1 = l_val / f_val * (cos(phi)*cos(phi)*g_val + sin(phi)*sin(phi));
    double g_xy_1 = (cos(phi)/f_val)*(-1+g_val)*l_val*sin(phi);
    // double g_yx_1 = (cos(phi)/f_val)*(-1+g_val)*l_val*sin(phi);

    //Add on to evolution equations
    vars.phi_Re += A_val * cos(phase_);
    vars.phi_Im += A_val * sin(phase_);
    vars.Pi_Re += (A_val / lapse) * (ff - beta_phi)*sin(phase_);
    vars.Pi_Im += -(A_val / lapse) * (ff - beta_phi)*cos(phase_);

    // pout() << "Computed real and imaginary scalar part expressions " << endl;

    //Initialise extrinsic curvature and metric with upper indices
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K;

    // Fill them in
    gammaLL[0][0] = g_xx_1;
    gammaLL[1][1] = g_yy_1;
    gammaLL[2][2] = g_zz_1;
    gammaLL[0][1] = g_xy_1;
    gammaLL[1][0] = gammaLL[0][1];

    gammaUU[0][0] = f_val / (g_val * l_val) * (cos(phi)*cos(phi) + g_val * sin(phi)*sin(phi));
    gammaUU[0][1] = (-cos(phi) * f_val * (-1. + g_val) * sin(phi)) / (g_val * l_val);
    gammaUU[1][1] = f_val / (g_val * l_val) * (cos(phi)*cos(phi) * g_val + sin(phi)*sin(phi));;
    gammaUU[1][0] = gammaUU[0][1];
    gammaUU[2][2] = 1. / g_zz_1;

    ///////////////////////
    //  For debugging    //
    ///////////////////////
    // double check11 = gammaLL[0][0]*gammaUU[0][0] + gammaLL[0][1]*gammaUU[1][0] + gammaLL[0][2] * gammaUU[2][0];
    // double check01 = gammaLL[0][0]*gammaUU[0][1] + gammaLL[0][1]*gammaUU[1][1] + gammaLL[0][2] * gammaUU[2][1];
    // double check02 = gammaLL[0][0]*gammaUU[0][2] + gammaLL[0][1]*gammaUU[1][2] + gammaLL[0][2] * gammaUU[2][2];
    // if (fabs(check11-1.0)>1e-3)
    // {   
    //     pout() << "check11 is not 1 but " << check11 << endl;
    // }
    // if (fabs(check01-0.0)>1e-3)
    // {   
    //     pout() << "check01 is not 0 but " << check01 << endl;
    // }
    // if (fabs(check02-0.0)>1e-3)
    // {   
    //     pout() << "check02 is not 0 but " << check02 << endl;
    // }

    KLL[0][0] = (l_val*sin(2.*phi)*(cos(theta) * sin(theta) * dthomega_val+sin(theta)*sin(theta)*(-omega_val + r*dromega_val)))/(4*f_val*r*lapse);
    KLL[0][1] = -(l_val*cos(2.*phi)*(cos(theta) * sin(theta) * dthomega_val+sin(theta)*sin(theta)*(-omega_val + r*dromega_val)))/(4*f_val*r*lapse);
    KLL[1][0] = KLL[0][1];
    KLL[1][1] = -KLL[0][0];
    KLL[0][2] = (-l_val*sin(phi)*(sin(theta)*sin(theta)*dthomega_val + cos(theta)*sin(theta)*(omega_val - r*dromega_val)))/(4.*f_val*r*lapse);
    KLL[2][0] = KLL[0][2];
    KLL[1][2] = (l_val*cos(phi)*(sin(theta)*sin(theta)*dthomega_val + cos(theta)*sin(theta)*(omega_val - r*dromega_val)))/(4.*f_val*r*lapse);
    KLL[2][1] = KLL[1][2];
    KLL[2][2] = 0;

    // pout() << "Computed Kij expressions " << endl;

    double chi_arg = (g_val * g_val * l_val * l_val * l_val  / (pow(f_val, 3)));
    vars.chi = pow(chi_arg, -1. / 3.);

    // pout() << "Computed BBSN conformal factor expressions " << endl;

    // Define initial lapse
    vars.lapse += lapse;

    // Define initial trace of K and A_ij
    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = vars.chi * gammaLL[i][j];
    // FOR2(i,j) 
    // {
    //     vars.K += KLL[i][j] * gammaUU[i][j];
    //     if (vars.K != 0)
    //     {
    //         pout() << "At i " << i << " and j " << j << "we have " << vars.K << endl;
    //     }
    // }
    // FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
    vars.K = 0.0;
    FOR2(i,j) vars.A[i][j] = vars.chi * KLL[i][j];

    current_cell.store_vars(vars);
}

#endif /* ROTATINGBOSONSTAR_IMPL_HPP_ */
