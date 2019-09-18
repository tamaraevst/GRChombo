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
    double separation = m_params_BosonStar.BS_separation;
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double t = coords.z*s_ + 0.*s_; //set /tilde{t} to zero
    double x = coords.x-separation/2.;
    double z = coords.z*c_;
    double y = coords.y;
    double r = sqrt(x*x+y*y+z*z);
    double r_tilde = sqrt((x+separation)*(x+separation)+y*y+z*z + 10e-10);
    double M = 5.327;
    double BH_lapse = (2.-M/r_tilde)/(2.+M/r_tilde);
    double BH_gxx = pow(1.+M/(2.*r_tilde),4);
    double p_ = m_1d_sol.get_p_interp(r);
    double dp_ = m_1d_sol.get_dp_interp(r);
    double alpha_ = m_1d_sol.get_lapse_interp(r);
    double alpha_prime_ = m_1d_sol.get_dlapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double psi_prime_ = m_1d_sol.get_dpsi_interp(r);
    double pc_os = psi_*psi_*c_*c_ - alpha_*alpha_*s_*s_;
    double lapse_ = alpha_*psi_/(sqrt(pc_os));
    double gamma_ = pow(psi_,4)*(pc_os);
    double w_ = m_1d_sol.get_w();
    double chi_ = pow(gamma_,-1./3.);
    double phase_ = m_params_BosonStar.phase + w_*t;
    double beta_z = s_*c_*(alpha_*alpha_-psi_*psi_)/(pc_os);
    double g_xx = psi_*psi_ + BH_gxx - 1.;
    double g_zz = pc_os + BH_gxx - 1.;
    double lapse_prime_ = (pow(psi_,3)*c_*c_*alpha_prime_ - pow(alpha_,3)*s_*s_*psi_prime_)/(pow(pc_os,1.5));
    double beta_z_prime_ = (psi_*alpha_*sinh(2.*rapidity)*(psi_*alpha_prime_-alpha_*psi_prime_))/(pc_os*pc_os);
    double L_n_gamma_1 = -(z/(lapse_*r))*(s_-beta_z*c_); // two coefficients of the lie derivative of 3 metric which nicely splits into 2 matrices
    double L_n_gamma_2 = (1./(2.*lapse_*lapse_*r))*(lapse_*beta_z_prime_-beta_z*lapse_prime_)*pc_os;
    vars.phi_Re += p_*cos(phase_);
    vars.phi_Im += p_*sin(phase_);
    vars.Pi_Re += -(1./lapse_)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
    vars.Pi_Im += -(1./lapse_)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );
    //vars.lapse += sqrt(lapse_*lapse_ + BH_lapse*BH_lapse - 1.);
    vars.shift[2] += beta_z;


    if (binary)
    {
        c_ = cosh(-rapidity);
        s_ = sinh(-rapidity);
        t = coords.z*s_ + 0.*s_; //set /tilde{t} to zero
        x = coords.x+separation/2.;
        z = coords.z*c_;
        y = coords.y;
        r = sqrt(x*x+y*y+z*z);
        p_ = m_1d_sol.get_p_interp(r);
        dp_ = m_1d_sol.get_dp_interp(r);
        alpha_ = m_1d_sol.get_lapse_interp(r);
        alpha_prime_ = m_1d_sol.get_dlapse_interp(r);
        psi_ = m_1d_sol.get_psi_interp(r);
        psi_prime_ = m_1d_sol.get_dpsi_interp(r);
        pc_os = psi_*psi_*c_*c_ - alpha_*alpha_*s_*s_;
        lapse_ = alpha_*psi_/(sqrt(pc_os));
        gamma_ = pow(psi_,4)*(pc_os);
        w_ = m_1d_sol.get_w();
        chi_ = pow(gamma_,-1./3.);
        phase_ = m_params_BosonStar.phase + w_*t;
        beta_z = s_*c_*(alpha_*alpha_-psi_*psi_)/(pc_os);
        g_xx += psi_*psi_-1.;
        g_zz += pc_os-1.;
        vars.phi_Re += p_*cos(phase_);
        vars.phi_Im += p_*sin(phase_);
        vars.Pi_Re += -(1./lapse_)*( (x/r)*(s_-beta_z*c_)*dp_*cos(phase_) - w_*(c_-beta_z*s_)*p_*sin(phase_) );
        vars.Pi_Im += -(1./lapse_)*( (x/r)*(s_-beta_z*c_)*dp_*sin(phase_) + w_*(c_-beta_z*s_)*p_*cos(phase_) );
        vars.lapse += lapse_-1.;
        vars.shift[2] += beta_z;
    }

    vars.chi += pow(g_xx*g_xx*g_zz,-1./3.);
    vars.lapse += sqrt(vars.chi);
    vars.h[0][0] += pow(g_xx*g_xx*g_zz,-1./3.)*g_xx;
    vars.h[1][1] += pow(g_xx*g_xx*g_zz,-1./3.)*g_xx;
    vars.h[2][2] += pow(g_xx*g_xx*g_zz,-1./3.)*g_zz;
    double K11 = L_n_gamma_1*psi_prime_*psi_;
    double K22 = L_n_gamma_1*psi_prime_*psi_;
    double K33 = L_n_gamma_1*(c_*c_*psi_*psi_prime_-s_*s_*lapse_*lapse_prime_) + L_n_gamma_2*2.*z*c_;
    double K13 = L_n_gamma_2*x ;
    double K23 = L_n_gamma_2*y ;
    double K12 = 0.;
    double one_third = 1./3.;
    vars.K += (K11+K22)/g_xx + K33/g_zz;
    vars.A[0][0] += chi_*(K11-one_third*vars.K*g_xx);
    vars.A[1][1] += chi_*(K22-one_third*vars.K*g_xx);
    vars.A[2][2] += chi_*(K33-one_third*vars.K*g_zz);
    vars.A[0][1] += chi_*K12;
    vars.A[0][2] += chi_*K13;
    vars.A[1][2] += chi_*K23;




    /*

    // define coords wrt first star centre
    // define coords wrt first star centre
    double rapidity = m_params_BosonStar.BS_rapidity;
    bool binary = m_params_BosonStar.BS_binary;
    double separation = m_params_BosonStar.BS_separation;
    double t = coords.z*sinh(rapidity);
    double x = coords.x-separation/2.;
    double z = coords.z*cosh(rapidity); //boosting star along y direction
    double y = coords.y;

    double r = sqrt(x*x+y*y+z*z);
    double p_ = m_1d_sol.get_p_interp(r);
    double dp_ = m_1d_sol.get_dp_interp(r);
    double alpha_ = m_1d_sol.get_lapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double lapse_ = psi_*alpha_/sqrt(pow(psi_*cosh(rapidity),2)-pow(alpha_*sinh(rapidity),2));
    double gamma_ = pow(psi_,4)*(pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2));
    double w_ = m_1d_sol.get_w();
    double chi_ = pow(gamma_,-1./3.);
    double phase_ = m_params_BosonStar.phase + w_*t;
    double g_yy = psi_*psi_;
    double g_zz = pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2);

    //Complex scalar field values
    vars.phi_Re += p_*cos(phase_);
    vars.phi_Im += p_*sin(phase_);
    vars.Pi_Re += -(psi_*cosh(rapidity)/(alpha_*sqrt(gamma_)))*( -w_*psi_*psi_*p_*sin(phase_)  +  ((coords.y/sqrt(r*r + 0.00001)))*alpha_*alpha_*dp_*sinh(rapidity)*cos(phase_) );
    vars.Pi_Im += -(psi_*cosh(rapidity)/(alpha_*sqrt(gamma_)))*( w_*psi_*psi_*p_*cos(phase_)  +  ((coords.y/sqrt(r*r + 0.00001)))*alpha_*alpha_*dp_*sinh(rapidity)*sin(phase_) );

    //conformal factor and lapse
    vars.chi += chi_;
    vars.lapse += lapse_;
    vars.shift[1] += sinh(rapidity)*cosh(rapidity)*(alpha_*alpha_-psi_*psi_)/( (pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2)) ); // beta^y shift

    if (binary)
    {
        // now superpose the second star
        t = coords.z*sinh(-rapidity); //
        x = coords.x+separation/2.;
        z = coords.z*cosh(-rapidity); //boosting star along -z
        y = coords.y;
        r = sqrt(x*x+y*y+z*z);
        p_ = m_1d_sol.get_p_interp(r);
        dp_ = m_1d_sol.get_dp_interp(r);
        alpha_ = m_1d_sol.get_lapse_interp(r);
        psi_ = m_1d_sol.get_psi_interp(r);
        lapse_ = psi_*alpha_/sqrt(pow(psi_*cosh(rapidity),2)-pow(alpha_*sinh(rapidity),2));
        gamma_ = pow(psi_,4)*(pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2));
        w_ = m_1d_sol.get_w(); // can make negative for opposite phase rotation
        chi_ = pow(gamma_,-1./3.);
        phase_ = m_params_BosonStar.phase + w_*t;
        g_zz += pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2)-1.;
        g_yy += psi_*psi_-1.;

        //Complex scalar field values
        vars.phi_Re += p_*cos(phase_);
        vars.phi_Im += p_*sin(phase_);
        vars.Pi_Re += -(psi_*cosh(-rapidity)/(alpha_*sqrt(gamma_)))*( -w_*psi_*psi_*p_*sin(phase_)  +  ((coords.y/sqrt(r*r + 0.00001)))*alpha_*alpha_*dp_*sinh(-rapidity)*cos(phase_) );
        vars.Pi_Im += -(psi_*cosh(-rapidity)/(alpha_*sqrt(gamma_)))*( w_*psi_*psi_*p_*cos(phase_)  +  ((coords.y/sqrt(r*r + 0.00001)))*alpha_*alpha_*dp_*sinh(-rapidity)*sin(phase_) );


        //conformal factor and lapse
        vars.chi += chi_-1;
        vars.lapse += lapse_-1;
        vars.shift[1] += sinh(rapidity)*cosh(rapidity)*(alpha_*alpha_-psi_*psi_)/( (pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2)) ); // beta^y shift
    }

    //conformal metric is flat
    chi_ = pow(g_yy*g_yy*g_zz,-1/3);
    vars.h[0][0] += chi_*g_yy; //g_yy = g_xx as we boost along z
    vars.h[1][1] += chi_*g_yy;
    vars.h[2][2] += chi_*g_zz;

    */

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
