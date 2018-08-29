/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARINTEGRATOR_HPP_)
#error "This file should only be included through BosonStarIntegrator.hpp"
#endif

#ifndef BOSONSTARINTEGRATOR_IMPL_HPP_
#define BOSONSTARINTEGRATOR_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarIntegrator<initial_data_t, initial_state_t>::BosonStarIntegrator(
    BosonStar::params_t a_params_BosonStar,
    Potential::params_t a_params_potential)
    : m_params_BosonStar(a_params_BosonStar),
    m_boson_star_rhs(a_params_potential),
    m_initial_var_arrays{},
    m_initial_grid{},
    m_sol_observer(m_initial_var_arrays, m_initial_grid) {}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIntegrator<initial_data_t, initial_state_t>::clearArrays()
{
    m_initial_var_arrays.clear();
    m_initial_grid.clear();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIntegrator<initial_data_t, initial_state_t>
    ::doIntegration(const double a_alpha_central)
{
    //First clear arrays to make sure nothing has been left from the last
    //integration.
    clearArrays();

    //identify fixed BCs
    const double beta_central{0.0};
    const double Psi_central{0.0};
    const double central_radius{0.0}; /*bit of a weird name but makes it easier
                                        to understand later on */
    //Set central BCs
    initial_state_t central_vars{a_alpha_central, beta_central,
        m_params_BosonStar.central_amplitude_CSF, Psi_central};

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<initial_state_t> error_stepper_t;
    //Need to put this in a try block as the solution observer can throw
    //exceptions.
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>
            (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
            m_boson_star_rhs, central_vars, central_radius,
            m_params_BosonStar.max_radius, m_params_BosonStar.initial_step_size,
            m_sol_observer);
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << " max radius = " <<
            m_initial_grid.back() << "\n";
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>
    BosonStarIntegrator<initial_data_t, initial_state_t>
    ::getSolution()
{
    return BosonStarSolution<initial_data_t, initial_state_t>
        (m_initial_var_arrays, m_initial_grid);
}


#endif /* BOSONSTARINTEGRATOR_IMPL_HPP_ */
