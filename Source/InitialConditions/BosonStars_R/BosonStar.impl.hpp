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

    // define coords wrt first star centre 
    double rapidity = m_params_BosonStar.BS_rapidity;
    bool binary = m_params_BosonStar.BS_binary;
    double separation = m_params_BosonStar.BS_separation;
    double t = coords.y*sinh(rapidity);
    double x = coords.x-separation/2.;
    double y = coords.y*cosh(rapidity); //boosting star along y direction
    double z = coords.z;
    
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
    double g_zz = psi_*psi_;
    double g_yy = pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2);

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
        t = coords.y*sinh(-rapidity); //
        x = coords.x+separation/2.;
        y = coords.y*cosh(-rapidity); //boosting star along -y
        z = coords.z;
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
        g_yy += pow(cosh(rapidity)*psi_,2)-pow(sinh(rapidity)*alpha_,2)-1.;
        g_zz += psi_*psi_-1.;

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
    chi_ = pow(g_zz*g_zz*g_yy,-1/3);
    vars.h[0][0] += chi_*g_zz; //g_zz = g_zz as we boost along y
    vars.h[1][1] += chi_*g_yy;
    vars.h[2][2] += chi_*g_zz;


    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
