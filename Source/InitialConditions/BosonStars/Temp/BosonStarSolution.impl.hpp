/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARSOLUTION_HPP_)
#error "This file should only be included through BosonStarSolution.hpp"
#endif

#ifndef BOSONSTARSOLUTION_IMPL_HPP_
#define BOSONSTARSOLUTION_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>::BosonStarSolution(
    initial_data_t<initial_state_t> &a_initial_var_arrays,
    initial_data_t<double> &a_radii) : m_radii(a_radii), m_alpha_array {},
    m_beta_array {}, m_psi_array {}, m_Psi_array {}, m_num_psi_roots(0)
    {
        separateArrays(a_initial_var_arrays);
        count_num_psi_roots();
    }
template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>::BosonStarSolution(
    initial_data_t<initial_state_t> &a_initial_var_arrays,
    initial_data_t<double> &a_radii, int a_num_psi_roots) : m_radii(a_radii),
    m_alpha_array {}, m_beta_array {}, m_psi_array {}, m_Psi_array {},
    m_num_psi_roots(a_num_psi_roots)
    {
        separateArrays(a_initial_var_arrays);
    }

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::separateArrays(
    initial_data_t<initial_state_t> &a_initial_var_arrays)
{
    //first get the number of grid values
     m_num_grid_values = m_radii.size();

    //allocate enough memory for the variable arrays
    m_alpha_array.reserve(m_num_grid_values);
    m_beta_array.reserve(m_num_grid_values);
    m_psi_array.reserve(m_num_grid_values);
    m_Psi_array.reserve(m_num_grid_values);

    //now put the values of m_initial_var_arrays into the separate variable
    //arrays
    for(int i = 0; i < m_num_grid_values; ++i)
    {
        m_alpha_array.push_back(a_initial_var_arrays[i][0]);
        m_beta_array.push_back(a_initial_var_arrays[i][1]);
        m_psi_array.push_back(a_initial_var_arrays[i][2]);
        m_Psi_array.push_back(a_initial_var_arrays[i][3]);
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>::
    get_max_radius() const
{
    return m_radii.back();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>::
    get_psi_at_max_radius() const
{
    return m_psi_array.back();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::count_num_psi_roots()
{
    for (int i = 1; i < m_num_grid_values; ++i)
    {
        if(m_psi_array[i - 1] * m_psi_array[i] < 0)
        {
            ++m_num_psi_roots;
        }
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
int BosonStarSolution<initial_data_t, initial_state_t>::get_num_psi_roots() const
{
    return m_num_psi_roots;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_grid()
{
    return m_radii;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_alpha()
{
    return m_alpha_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_beta()
{
    return m_beta_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_psi()
{
    return m_psi_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_Psi()
{
    return m_Psi_array;
}

#endif /* BOSONSTARSOLUTION_IMPL_HPP_ */
