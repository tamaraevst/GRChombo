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
    double t = coords.z*s_ + 0.*s_; //set /tilde{t} to zero
    double x = coords.x-separation/2.;
    double z = coords.z*c_;
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
    double chi_;
    double phase_ = m_params_BosonStar.phase + w_*t;
    double beta_z = s_*c_*(psi_*psi_-omega_*omega_)/(pc_os);
    double g_xx = psi_*psi_;
    double g_yy = psi_*psi_;
    double g_zz = pc_os;

    vars.phi_Re += p_*cos(phase_);
    vars.phi_Im += p_*sin(phase_);
    vars.Pi_Re += -(1./lapse_1)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
    vars.Pi_Im += -(1./lapse_1)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );

    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double kLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K;

    kLL[0][0] = -lapse_1*s_*z*psi_prime_/(r*psi_);
    KLL[1][1] = KLL[0][0];
    KLL[2][1] = lapse_1*c_*s_*(y/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL[2][0] = lapse_1*c_*s_*(x/r)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL[1][2] = KLL[2][1];
    KLL[0][2] = KLL[2][0];
    KLL[0][1] = 0.;
    KLL[1][0] = 0.;
    KLL[2][2] = lapse_1*(z/r)*s_*c_*c_*(psi_prime_/psi_ - 2.*omega_prime_/omega_ + v_*v_*omega_*omega_prime_*pow(psi_,-2));
    double thingy = omega_*psi_prime_*(2.*v_*v_*omega_*omega_-psi_*psi_) + omega_prime_*psi_*(v_*v_*omega_*omega_-2.*psi_*psi_);
    K = pow(lapse_1/omega_,3)*pow(psi_,-5)*s_*c_*c_*(x/r)*(thingy);

    if (binary)
    {
        c_ = cosh(rapidity);
        s_ = sinh(rapidity);
        v_ = tanh(rapidity);
        t = coords.z*s_ + 0.*s_; //set /tilde{t} to zero
        x = coords.x-separation/2.;
        z = coords.z*c_;
        y = coords.y;
        r = sqrt(x*x+y*y+z*z);

          // first star physical variables
        p_ = m_1d_sol.get_p_interp(r);
        dp_ = m_1d_sol.get_dp_interp(r);
        omega_ = m_1d_sol.get_lapse_interp(r);
        omega_prime_ = m_1d_sol.get_dlapse_interp(r);
        psi_ = m_1d_sol.get_psi_interp(r);
        psi_prime_ = m_1d_sol.get_dpsi_interp(r);

        pc_os = psi_*psi_*c_*c_ - omega_*omega_*s_*s_;
        lapse_2 = omega_*psi_/(sqrt(pc_os));
        w_ = m_1d_sol.get_w();
        phase_ = m_params_BosonStar.phase + w_*t;
        beta_z = s_*c_*(psi_*psi_-omega_*omega_)/(pc_os);
        g_xx += psi_*psi_-1.;
        g_yy += psi_*psi_-1.;
        g_zz += pc_os-1.;

        vars.phi_Re += p_*cos(phase_);
        vars.phi_Im += p_*sin(phase_);
        vars.Pi_Re += -(1./lapse_2)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
        vars.Pi_Im += -(1./lapse_2)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );

    }
    else if (BS_BH_binary)
    {
        double r_tilde = sqrt((x+separation)*(x+separation)+y*y+z*z + 10e-10);
        double BH_lapse = (2.-M/r_tilde)/(2.+M/r_tilde);
        double BH_gii = pow(1.+M/(2.*r_tilde),4);
        g_xx += BH_gii-1.;
        g_yy += BH_gii-1.;
        g_zz += BH_gii-1.;
    }

    vars.chi += pow(g_xx*g_yy*g_zz,-1./3.);
    vars.lapse += BS_BH_binary?sqrt(vars.chi):sqrt(lapse_1*lapse_1 + lapse_2+lapse_2-1.);
    vars.h[0][0] += vars.chi*g_xx;
    vars.h[1][1] += vars.chi*g_yy;
    vars.h[2][2] += vars.chi*g_zz;


    double one_third = 1./3.;
    vars.K += K;
    vars.A[0][0] += chi_*(KLL[0][0]-one_third*K*g_xx);
    vars.A[1][1] += chi_*(KLL[1][1]-one_third*K*g_yy);
    vars.A[2][2] += chi_*(KLL[2][2]-one_third*K*g_zz);
    vars.A[0][1] += chi_*KLL[0][1];
    vars.A[0][2] += chi_*KLL[0][2];
    vars.A[1][2] += chi_*KLL[1][2];
    vars.A[1][0] += chi_*KLL[1][0];
    vars.A[2][0] += chi_*KLL[2][0];
    vars.A[2][1] += chi_*KLL[2][1];

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
