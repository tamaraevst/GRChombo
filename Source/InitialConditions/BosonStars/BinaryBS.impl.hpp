/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBS_HPP_)
#error "This file should only be included through BinaryBS.impl.hpp"
#endif

#ifndef BINARYBS_IMPL_HPP_
#define BINARYBS_IMPL_HPP_

inline BinaryBS::BinaryBS(BosonStar_params_t a_boson_star1_params,
                          BosonStar_params_t a_boson_star2_params,
                          Potential::params_t a_params_potential,
                          double a_G_Newton, double a_dx, bool a_identical,
                          int a_verbosity)
    : m_dx(a_dx), m_boson_star1(a_boson_star1_params, a_params_potential,
        a_G_Newton, a_dx, a_verbosity), m_boson_star2(a_boson_star2_params,
        a_params_potential, a_G_Newton, a_dx, a_verbosity),
        m_identical(a_identical), m_verbosity(a_verbosity)
{
}

void BinaryBS::compute_profiles(const double a_L)
{
    if(m_verbosity)
    {
        pout() << "BinaryBS::compute_profiles: Computing boson star 1 profile"
               << std::endl;
    }
    m_boson_star1.compute_1d_solution(3.5 * a_L);

    // only need to compute profile for 2nd star if different from first
    if(!m_identical)
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles:"
                      " Computing boson star 2 profile" << std::endl;
        }
        m_boson_star2.compute_1d_solution(3.5 * a_L);
    }
    else
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles: Boson star 2 identical"
                      " to star 1; skipping profile computation" << std::endl;
        }

        // Copy boson_star1 into boson_star2 keeping phase and centre
        auto boson_star2_centre = m_boson_star2.m_params_BosonStar.star_centre;
        double boson_star2_phase = m_boson_star2.m_params_BosonStar.phase;
        m_boson_star2 = m_boson_star1;
        m_boson_star2.m_params_BosonStar.star_centre = boson_star2_centre;
        m_boson_star2.m_params_BosonStar.phase = boson_star2_phase;
    }
}

template <class data_t>
void BinaryBS::compute(Cell<data_t> current_cell) const
{
    // First just superpose all variables
    m_boson_star1.compute(current_cell);
    m_boson_star2.compute(current_cell);

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
