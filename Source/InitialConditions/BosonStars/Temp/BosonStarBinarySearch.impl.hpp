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
    BosonStar::params_t a_params_CSF, Potential::params_t a_params_potential,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_min,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_max)
    : m_params_CSF(a_params_CSF), m_params_potential(a_params_potential),
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
void BosonStarBinarySearch<initial_data_t, initial_state_t>
    ::shoot(const double a_tol, const double a_max_radius)
{
    //identify the two BCs that are fixed
    double beta_central{0.0};
    double Psi_central{0.0};

    //Set integration error tolerances
    double abs_error{1.0e-14};
    double rel_error{1.0e-14};

    //Set initial step size
    double drho{2.0e-7};

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
    typedef runge_kutta_cash_karp54<initial_state_t> error_stepper_t;

    //keep bisecting until we are within the desired tolerance
    while(m_alpha_central_max - m_alpha_central_min > a_tol)
    {
        //Clear the storage arrays from the last iteration
        initial_var_arrays_mid.clear();
        initial_grid_mid.clear();

        //Set central BCs for this iteration
        m_alpha_central_mid =
            0.5 * ( m_alpha_central_min + m_alpha_central_max );
        initial_state_t central_vars_mid{m_alpha_central_mid, beta_central,
            m_params_CSF.central_amplitude_CSF, Psi_central};

        //do integration
        try
        {
            integrate_adaptive(make_controlled<error_stepper_t>(
                abs_error, rel_error), boson_star_rhs, central_vars_mid, 0.0,
                a_max_radius, drho, sol_observer);
        }
        catch (std::exception &exception)
        {
            //std::cout << exception.what() << " max radius = " <<
            //    initial_grid_mid.back() << "\n";
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
        m_params_CSF.central_amplitude_CSF, Psi_central};

    //do integration for final mid solution
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>(
            abs_error, rel_error), boson_star_rhs, central_vars_mid, 0.0,
            a_max_radius, drho, sol_observer);
    }
    catch (std::exception &exception)
    {
        //std::cout << exception.what() << " max radius = " <<
        //    initial_grid_mid.back() << "\n";
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
