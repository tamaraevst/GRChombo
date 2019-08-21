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
    double x = coords.x-32;
    double y = coords.y;
    double z = coords.z;
    
    double r = sqrt(x*x+y*y+z*z);
    double p_ = m_1d_sol.get_p_interp(r);
    double lapse_ = m_1d_sol.get_lapse_interp(r);
    double w_ = m_1d_sol.get_w();
    double chi_ = m_1d_sol.get_chi_interp(r);
    double phi_ = m_params_BosonStar.phase;

    //Complex scalar field values
    vars.phi_Re += p_*cos(phi_);
    vars.phi_Im += p_*sin(phi_);
    vars.Pi_Re += p_*sin(phi_)*w_/lapse_;
    vars.Pi_Im += -p_*cos(phi_)*w_/lapse_;

    //conformal factor and lapse
    vars.chi += chi_;
    vars.lapse += lapse_;

    // now superpose the second star
    x += 64;
    r = sqrt(x*x+y*y+z*z);
    p_ = m_1d_sol.get_p_interp(r);
    lapse_ = m_1d_sol.get_lapse_interp(r);
    w_ = m_1d_sol.get_w(); // can make w negative for opposite phase rotation
    chi_ = m_1d_sol.get_chi_interp(r);
    phi_ = m_params_BosonStar.phase;

    //Complex scalar field values
    vars.phi_Re += p_*cos(phi_);
    vars.phi_Im += p_*sin(phi_);
    vars.Pi_Re += p_*sin(phi_)*w_/lapse_;
    vars.Pi_Im += -p_*cos(phi_)*w_/lapse_;
    
    //conformal factor and lapse
    vars.chi += chi_;
    vars.lapse += lapse_;

    //conformal metric is flat
    FOR1(i) vars.h[i][i] += 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
