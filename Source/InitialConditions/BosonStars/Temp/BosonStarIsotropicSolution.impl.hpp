/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARISOTROPICSOLUTION_HPP_)
#error "This file should only be included through BosonStarIsotropicSolution.hpp"
#endif

#ifndef BOSONSTARISOTROPICSOLUTION_IMPL_HPP_
#define BOSONSTARISOTROPICSOLUTION_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::BosonStarIsotropicSolution(
    BosonStarSolution<initial_data_t, initial_state_t>
    &a_polar_areal_solution, BosonStar::params_t a_params_BosonStar,
    Potential::params_t a_params_potential,const double a_max_radius)
    : m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential),
    m_polar_areal_solution(a_polar_areal_solution)
{
    calculateIsotropicGrid(a_max_radius);
    std::cout.precision(16);
    std::cout << std::fixed;
    std::cout << "rho\t\t\tR\n";
    for(int i = 0; i < static_cast<int>(m_polar_areal_grid.size()); ++i)
    {
        std::cout << m_polar_areal_grid[i] << "\t" << m_isotropic_grid[i]
        << "\n";
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::calculateIsotropicGrid(const double a_max_radius)
{
    //First interpolate beta
    tools::spline<initial_data_t> beta_interpolated;
    beta_interpolated.set_points(m_polar_areal_solution.get_grid(),
        m_polar_areal_solution.get_beta());

    //RHS for converting to isotropic coordinates
    auto rhs_lambda = [&](const initial_state_t &R , initial_state_t &dRdr ,
        double r )
    {
        dRdr[0] = std::exp(beta_interpolated(m_params_potential.scalar_mass \
            * r)) * R[0]/r;
    };

    //outer boundary condition
    initial_state_t R_max {0.25 * a_max_radius *
    (std::exp(-beta_interpolated(m_params_potential.scalar_mass * \
        a_max_radius)) + 1.0) * (std::exp(-beta_interpolated(\
            m_params_potential.scalar_mass * a_max_radius)) + 1.0)};

    //Since the ODE involves division by r, we can only integrate up to some
    //r > 0
    const double r_inner_limit {1.0e-10};

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<initial_state_t> error_stepper_t;
    integrate_adaptive(make_controlled<error_stepper_t>
        (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
        rhs_lambda, R_max, a_max_radius, r_inner_limit,
        -m_params_BosonStar.initial_step_size, //need to reverse sign as integrating inwards
        [&] (const initial_state_t &a_R, double &a_r)
        {
            m_isotropic_grid.push_back(a_R[0]);
            m_polar_areal_grid.push_back(a_r);
        });

    //since we can't integrate to the origin, add it in manually
    m_isotropic_grid.push_back(0.0);
    m_polar_areal_grid.push_back(0.0);

    //reverse both arrays as they currently start from the outer limit
    std::reverse(m_isotropic_grid.begin(), m_isotropic_grid.end());
    std::reverse(m_polar_areal_grid.begin(), m_polar_areal_grid.end());
}

#endif /* BOSONSTARISOTROPICSOLUTION_IMPL_HPP_ */
