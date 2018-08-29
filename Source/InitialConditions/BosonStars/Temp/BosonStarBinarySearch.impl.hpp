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
    : BosonStarIntegrator<initial_data_t, initial_state_t>
    (a_params_BosonStar, a_params_potential),
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

//Note that for some reason, C++ won't let you access member variables of a
//template base class without derefrencing the implicit this pointer. If this
//needs to be modified and any member function/variable of BosonStarIntegrator
//need to be accessed, remember to derefrence them from this*.
template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarBinarySearch<initial_data_t, initial_state_t>::shoot()
{
    int n_binary_search_iter{0};

    //keep bisecting until we are within the desired tolerance
    while(m_alpha_central_max - m_alpha_central_min
        > this->m_params_BosonStar.binary_search_tol &&
        n_binary_search_iter <= this->m_params_BosonStar.max_binary_search_iter)
    {
        //Set central BCs for this iteration
        m_alpha_central_mid =
            0.5 * ( m_alpha_central_min + m_alpha_central_max );
        this->doIntegration(m_alpha_central_mid);

        //Copy the storage arrays into a solution object
        m_sol_mid = this->getSolution();

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

    //Calculate final mid solution.
    m_alpha_central_mid = 0.5 * ( m_alpha_central_min + m_alpha_central_max );
    this->doIntegration(m_alpha_central_mid);

    //Copy the storage arrays into a solution object
    m_sol_mid = this->getSolution();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>&
    BosonStarBinarySearch<initial_data_t, initial_state_t>::getShootedSolution()
{
    return m_sol_mid;
}

#endif /* BOSONSTARBINARYSEARCH_IMPL_HPP_ */
