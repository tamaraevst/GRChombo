/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARBINARYSEARCH_HPP_)
#error "This file should only be included through BosonStarBinarySearch.hpp"
#endif

#ifndef BOSONSTARBINARYSEARCH_IMPL_HPP_
#define BOSONSTARBINARYSEARCH_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarBinarySearch<initial_data_t, initial_state_t>::BosonStarBinarySearch(
    BosonStar::params_t a_params_BosonStar,
    Potential::params_t a_params_potential,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_min,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_max)
    : m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential),
    m_alpha_central_min(a_sol_min.get_alpha()[0]),
    m_alpha_central_max(a_sol_max.get_alpha()[0]),
    m_sol_min(a_sol_min), m_sol_max(a_sol_max), m_sol_mid(a_sol_min)
{
    if( !checkValidInterval() )
    {
        throw std::runtime_error("The interval provided does not bound a"\
        " ground state soluion. Please try again!");
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
bool BosonStarBinarySearch<initial_data_t, initial_state_t>
    ::checkValidInterval()
{
    return( m_sol_min.get_psi_at_max_radius() * m_sol_max.get_psi_at_max_radius()
        < 0.0 && m_sol_min.get_num_psi_roots() == 1
        && m_sol_max.get_num_psi_roots() == 0 );
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarBinarySearch<initial_data_t, initial_state_t>::shoot()
{
    //identify the two BCs that are fixed
    double beta_central{0.0};
    double Psi_central{0.0};

    //initialise storage arrays
    initial_data_t<initial_state_t> initial_var_arrays_mid{};
    initial_grid_t initial_grid_mid{};

    //initialise RHS Class and solution observer. Since odeint makes a copy of
    //sol_observer, it should be safe to just instantiate it once here for all
    //the iterations in the while loop.
    BosonStarRHS boson_star_rhs(m_params_potential);
    BosonStarSolutionObserver<initial_data_t, initial_state_t>
        sol_observer(initial_var_arrays_mid, initial_grid_mid);

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<initial_state_t> error_stepper_t;

    //count the number of iterations so the loop doesn't run away
    int n_binary_search_iter{0};

    //keep bisecting until we are within the desired tolerance
    while(m_alpha_central_max - m_alpha_central_min
        > m_params_BosonStar.binary_search_tol &&
        n_binary_search_iter <= m_params_BosonStar.max_binary_search_iter)
    {
        //Clear the storage arrays from the last iteration
        initial_var_arrays_mid.clear();
        initial_grid_mid.clear();

        //Set central BCs for this iteration
        m_alpha_central_mid =
            0.5 * ( m_alpha_central_min + m_alpha_central_max );
        initial_state_t central_vars_mid{m_alpha_central_mid, beta_central,
            m_params_BosonStar.central_amplitude_CSF, Psi_central};

        //do integration
        try
        {
            integrate_adaptive(make_controlled<error_stepper_t>
                (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
                boson_star_rhs, central_vars_mid, 0.0,
                m_params_BosonStar.max_radius,
                m_params_BosonStar.initial_step_size, sol_observer);
        }
        catch (std::exception &exception)
        {
            std::cout << exception.what() << " max radius = " <<
                initial_grid_mid.back() << "\n";
        }

        //Copy the storage arrays into a solution object
        m_sol_mid = BosonStarSolution<initial_data_t, initial_state_t>
            (initial_var_arrays_mid, initial_grid_mid);

        if(m_sol_min.get_psi_at_max_radius()
            * m_sol_mid.get_psi_at_max_radius() > 0.0)
        {
            m_sol_min = m_sol_mid;
            m_alpha_central_min = m_alpha_central_mid;
        }
        else
        {
            m_sol_max = m_sol_mid;
            m_alpha_central_max = m_alpha_central_mid;
        }
    }

    //clear storage arrays from last iteration
    initial_var_arrays_mid.clear();
    initial_grid_mid.clear();

    //Calculate final mid solution.
    m_alpha_central_mid = 0.5 * ( m_alpha_central_min + m_alpha_central_max );
    initial_state_t central_vars_mid{m_alpha_central_mid, beta_central,
        m_params_BosonStar.central_amplitude_CSF, Psi_central};

    //do integration for final mid solution
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>
            (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
            boson_star_rhs, central_vars_mid, 0.0,
            m_params_BosonStar.max_radius,
            m_params_BosonStar.initial_step_size, sol_observer);
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << " max radius = " <<
            initial_grid_mid.back() << "\n";
    }
    //Copy the storage arrays into a solution object
    m_sol_mid = BosonStarSolution<initial_data_t, initial_state_t>
        (initial_var_arrays_mid, initial_grid_mid);
}

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>&
    BosonStarBinarySearch<initial_data_t, initial_state_t>::getSolution()
{
    return m_sol_mid;
}

#endif /* BOSONSTARBINARYSEARCH_IMPL_HPP_ */
