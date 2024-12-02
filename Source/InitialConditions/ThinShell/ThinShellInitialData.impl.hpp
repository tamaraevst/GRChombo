/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(THINSHELLINITIALDATA_HPP_)
#error "This file should only be included through ThinShellInitialData.hpp"
#endif

#ifndef THINSHELLINITIALDATA_IMPL_HPP_
#define THINSHELLINITIALDATA_IMPL_HPP_

#include <cmath>
#include "ThinShellInitialData.hpp"
#include "DebuggingTools.hpp"

inline ThinShellInitialData::ThinShellInitialData(ThinShell_params_t a_params_ThinShell, double a_dx)
    : m_params_ThinShell(a_params_ThinShell), m_dx(a_dx)
{
}

void ThinShellInitialData::compute_1d_solution()
{   
    try
    {  
        thin_shell_sol.main(m_params_ThinShell.base_path);
        pout() << "Wooo I have read the spinning BS initial data!" << endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ThinShellInitialData::compute(Cell<data_t> current_cell) const
{   
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_ThinShell.star_centre);

    // Star positioning
    double t = 0;
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);

    double ff = m_params_ThinShell.BS_frequency;

    int m = thin_shell_sol.m;
    
    double A_val = thin_shell_sol.get_amp_interp(r);
    double f_val = thin_shell_sol.get_f_interp(r);
    double Phi_val = thin_shell_sol.get_lapse_interp(r);

    double lapse = exp(Phi_val);
    double beta_x = 0;
    double beta_y = 0;
    double beta_z = 0;

    vars.shift[0] += beta_x;
    vars.shift[1] += beta_y;
    vars.shift[2] += beta_z;

    // pout() << "Computed lapse and shift expressions \n" << endl;

    double phase_ = ff * t;

    double g_xx_1 = pow(f_val, -2);
    double g_yy_1 = pow(f_val, -2);
    double g_zz_1 = pow(f_val, -2);

    //Add on to evolution equations
    vars.phi_Re += A_val;
    vars.phi_Im += 0.0;
    vars.Pi_Re += 0.0;
    vars.Pi_Im += -(1. /lapse) * (ff * A_val);

    //Initialise extrinsic curvature and metric with upper indices
    double KLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K1;

    gammaLL_1[0][0] = g_xx_1;
    gammaLL_1[1][1] = g_yy_1;
    gammaLL_1[2][2] = g_zz_1;

    vars.chi = pow(g_xx_1 * g_yy_1 * g_zz_1, -1. / 3.);

    // Define initial lapse
    vars.lapse += lapse;

    // Define initial trace of K and A_ij
    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = vars.chi * gammaLL_1[i][j];

    vars.K = 0.0;

    FOR2(i,j) vars.A[i][j] = 0.0;

    current_cell.store_vars(vars);
}

#endif /* ROTATINGBOSONSTAR_IMPL_HPP_ */
