/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_IMPL_HPP_
#define BOSONSTAR_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class

inline BosonStar::BosonStar(BosonStar_params_t a_params_BosonStar,
                    Potential::params_t a_params_potential, double a_G_Newton,
                    double a_dx, int a_verbosity)
    :m_dx(a_dx), m_G_Newton(a_G_Newton), m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential), m_verbosity(a_verbosity)
{
}

void BosonStar::compute_1d_solution(const double max_r)
{
    try
    {
        //set initial parameters and then run the solver (didnt put it in the constructor)
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,m_params_potential,max_r);
        m_1d_sol.main();
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BosonStar::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);


    double rapidity = m_params_BosonStar.BS_rapidity;
    bool binary = m_params_BosonStar.BS_binary;
    bool BS_BH_binary = m_params_BosonStar.BS_BH_binary;
    double M = m_params_BosonStar.BlackHoleMass;
    double separation = m_params_BosonStar.BS_separation;

    // boosts and coordinate objects
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double v_ = tanh(rapidity);
    double t = coords.z*s_; //set /tilde{t} to zero
    double x = coords.x-separation/2.;
    double z = (coords.z+0.*separation/2.)*c_; //set /tilde{t} to zero
    double y = coords.y;
    double r = sqrt(x*x+y*y+z*z);

    // first star physical variables
    double p_ = m_1d_sol.get_p_interp(r);
    double dp_ = m_1d_sol.get_dp_interp(r);
    double omega_ = m_1d_sol.get_lapse_interp(r);
    double omega_prime_ = m_1d_sol.get_dlapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double psi_prime_ = m_1d_sol.get_dpsi_interp(r);

    double pc_os = psi_*psi_*c_*c_ - omega_*omega_*s_*s_;
    double lapse_1 = omega_*psi_/(sqrt(pc_os));
    double lapse_2 = 1.;
    double w_ = m_1d_sol.get_w();
    double phase_ = w_*t;
    double beta_z = s_*c_*(psi_*psi_-omega_*omega_)/(pc_os);
    vars.shift[2] += beta_z;
    double g_xx_1 = psi_*psi_;
    double g_yy_1 = psi_*psi_;
    double g_zz_1 = pc_os;
    double g_xx_2=1., g_yy_2=1., g_zz_2=1., g_xx, g_yy, g_zz;

    vars.phi_Re += p_*cos(phase_);
    vars.phi_Im += p_*sin(phase_);
    vars.Pi_Re += -(1./lapse_1)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
    vars.Pi_Im += -(1./lapse_1)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );

    double KLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K1, K2;
    gammaUU_1[0][0] = 1./g_xx_1;
    gammaUU_1[1][1] = 1./g_yy_1;
    gammaUU_1[2][2] = 1./g_zz_1;

    KLL_1[0][0] = -lapse_1*s_*z*psi_prime_/(r*psi_);
    KLL_1[1][1] = KLL_1[0][0];
    KLL_1[2][1] = lapse_1*c_*s_*(y/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL_1[2][0] = lapse_1*c_*s_*(x/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL_1[1][2] = KLL_1[2][1];
    KLL_1[0][2] = KLL_1[2][0];
    KLL_1[0][1] = 0.;
    KLL_1[1][0] = 0.;
    KLL_1[2][2] = lapse_1*(z/r)*s_*c_*c_*(psi_prime_/psi_ - 2.*omega_prime_/omega_ + v_*v_*omega_*omega_prime_*pow(psi_,-2));
    double thingy1 = omega_*psi_prime_*(2.*v_*v_*omega_*omega_-psi_*psi_) + omega_prime_*psi_*(v_*v_*omega_*omega_-2.*psi_*psi_);
    K1 = pow(lapse_1/omega_,3)*pow(psi_,-5)*s_*c_*c_*(z/r)*(thingy1);

    if (binary)
    {
        c_ = cosh(-rapidity);
        s_ = sinh(-rapidity);
        v_ = tanh(-rapidity);
        t = coords.z*s_; //set /tilde{t} to zero
        x = coords.x+separation/2.;
        z = (coords.z-0.*separation/2.)*c_;
        y = coords.y;
        r = sqrt(x*x+y*y+z*z);

          // first star physical variables
        p_ = m_1d_sol.get_p_interp(r);
        dp_ = m_1d_sol.get_dp_interp(r);
        omega_ = m_1d_sol.get_lapse_interp(r);
        omega_prime_ = m_1d_sol.get_dlapse_interp(r);
        psi_ = m_1d_sol.get_psi_interp(r);
        psi_prime_ = m_1d_sol.get_dpsi_interp(r);
        double r_tilde;

        if (BS_BH_binary)
        {
            r_tilde = sqrt(r*r + 10e-10);
            omega_ = (2.-M/r_tilde)/(2.+M/r_tilde);
            omega_prime_ = 4.*M/pow(2.*r_tilde + M,2);
            psi_ = pow(1.+M/(2.*r_tilde),2);
            psi_prime_ = -(M/(r_tilde*r_tilde))*(1.+M/(2.*r_tilde));
        }

        pc_os = psi_*psi_*c_*c_ - omega_*omega_*s_*s_;
        lapse_2 = omega_*psi_/(sqrt(pc_os));
        w_ = m_1d_sol.get_w();
        phase_ = m_params_BosonStar.phase*M_PI + w_*t;
        beta_z = s_*c_*(psi_*psi_-omega_*omega_)/(pc_os);
        vars.shift[2] += beta_z;
        g_xx_2 = psi_*psi_;
        g_yy_2 = psi_*psi_;
        g_zz_2 = pc_os;
        gammaUU_2[0][0] = 1./g_xx_2;
        gammaUU_2[1][1] = 1./g_yy_2;
        gammaUU_2[2][2] = 1./g_zz_2;


        if (BS_BH_binary)
        {
            //do not need to modify scalar field or momentum if we have a black hole
        }
        else
        {
            vars.phi_Re += p_*cos(phase_);
            vars.phi_Im += p_*sin(phase_);
            vars.Pi_Re += -(1./lapse_2)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
            vars.Pi_Im += -(1./lapse_2)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );
        }



        KLL_2[0][0] = -lapse_2*s_*z*psi_prime_/(r*psi_);
        KLL_2[1][1] = KLL_2[0][0];
        KLL_2[2][1] = lapse_2*c_*s_*(y/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
        KLL_2[2][0] = lapse_2*c_*s_*(x/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
        KLL_2[1][2] = KLL_2[2][1];
        KLL_2[0][2] = KLL_2[2][0];
        KLL_2[0][1] = 0.;
        KLL_2[1][0] = 0.;
        KLL_2[2][2] = lapse_2*(z/r)*s_*c_*c_*(psi_prime_/psi_ - 2.*omega_prime_/omega_ + v_*v_*omega_*omega_prime_*pow(psi_,-2));
        double thingy2 = omega_*psi_prime_*(2.*v_*v_*omega_*omega_-psi_*psi_) + omega_prime_*psi_*(v_*v_*omega_*omega_-2.*psi_*psi_);
        K2 = pow(lapse_2/omega_,3)*pow(psi_,-5)*s_*c_*c_*(z/r)*(thingy2);

    }
    g_xx = g_xx_1 + g_xx_2 - 1.;
    g_yy = g_yy_1 + g_yy_2 - 1.;
    g_zz = g_zz_1 + g_zz_2 - 1.;
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1./g_xx;
    gammaUU[1][1] = 1./g_yy;
    gammaUU[2][2] = 1./g_zz;

    double chi_ = pow(g_xx*g_yy*g_zz,-1./3.);
    vars.chi = chi_;


    if (BS_BH_binary){vars.lapse += sqrt(vars.chi);}
    else if (binary){vars.lapse += sqrt(lapse_1*lapse_1 + lapse_2*lapse_2-1.);}
    else{vars.lapse += lapse_1;}

    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = vars.chi*gammaLL[i][j];
    FOR2(i,j) vars.K += (KLL_1[i][j]+KLL_2[i][j])*gammaUU[i][j];
    FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l]*(gammaUU_1[l][k]*KLL_1[k][j] + gammaUU_2[l][k]*KLL_2[k][j]);
    FOR2(i,j) vars.A[i][j] = chi_*(KLL[i][j]-one_third*vars.K*gammaLL[i][j]);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
