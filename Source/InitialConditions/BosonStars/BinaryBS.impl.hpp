/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBS_HPP_)
#error "This file should only be included through BinaryBS.impl.hpp"
#endif

#ifndef BINARYBS_IMPL_HPP_
#define BINARYBS_IMPL_HPP_

inline BinaryBS::BinaryBS(BosonStar_params_t a_params_BosonStar1,
                          BosonStar_params_t a_params_BosonStar2,
                          Potential::params_t a_params_potential,
                          double a_G_Newton, double a_dx, bool a_identical,
                          int a_verbosity)
    : m_dx(a_dx), m_BosonStar1(a_params_BosonStar1, a_params_potential,
        a_G_Newton, a_verbosity), m_BosonStar2(a_params_BosonStar2,
        a_params_potential, a_G_Newton, a_verbosity), m_identical(a_identical),
        m_verbosity(a_verbosity)
{
}

void BinaryBS::compute_profiles(const double a_L)
{
    if(m_verbosity)
    {
        pout() << "BinaryBS::compute_profiles: Computing boson star 1 profile"
               << std::endl;
    }
    m_BosonStar1.compute_1d_solution(3.5 * a_L);

    // only need to compute profile for 2nd star if different from first
    if(!m_identical)
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles:"
                      " Computing boson star 2 profile" << std::endl;
        }
        m_BosonStar2.compute_1d_solution(3.5 * a_L);
    }
    else
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles: Boson star 2 identical"
                      " to star 1; skipping profile computation" << std::endl;
        }
    }
}

template <class data_t>
void BinaryBS::compute(Cell<data_t> current_cell) const
{
    // First just superpose all variables
    m_BosonStar1.compute(current_cell);
    if(!m_identical)
    {
        m_BosonStar2.compute(current_cell);
    }
    else
    {
        m_BosonStar1.m_params_BosonStar.star_centre
            = m_BosonStar2.m_params_BosonStar.star_centre;
        m_BosonStar1.m_params_BosonStar.phase
            = m_BosonStar2.m_params_BosonStar.phase;
        m_BosonStar1.compute(current_cell);
    }
    // Now manipulate here
    CCZ4Vars::VarsWithGauge<data_t> vars;
    current_cell.load_vars(vars);

    // Substract off Minkowski metric
    FOR1(i) vars.h[i][i] -= 1.;
    vars.chi -= 1.;
    vars.lapse -= 1.;

    // Store variables
    current_cell.store_vars(vars);
}

#endif /* BINARYBS_IMPL_HPP_ */
