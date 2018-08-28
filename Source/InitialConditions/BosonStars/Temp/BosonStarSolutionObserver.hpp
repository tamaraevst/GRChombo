/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARSOLUTIONOBSERVER_HPP_
#define BOSONSTARSOLUTIONOBSERVER_HPP_

#include <cmath> //for std::abs
#include <stdexcept> //for std::exception

//! This is a class to store the solution found by odeint in storage arrays as
//! it is calculated.

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarSolutionObserver
{
public:
    initial_data_t<initial_state_t> &m_initial_var_arrays; /*!< arrays that
        hold the grid values for all rescaled initial quantities */
    initial_data_t<double> &m_initial_grid; /*!< array that stores the rescaled radial
        coordinate corresponding to the grid values */
    int m_num_psi_roots; //!< this counts the number of roots in psi

    //! Constructor
    BosonStarSolutionObserver(
        initial_data_t<initial_state_t> &a_initial_var_arrays,
        initial_data_t<double> &a_initial_grid, int a_num_psi_roots = 0)
         : m_initial_var_arrays(a_initial_var_arrays),
         m_initial_grid(a_initial_grid), m_num_psi_roots(a_num_psi_roots) {}

    //! Overloaded () operator required for interface with boost odeint library
    void operator() (const initial_state_t &a_vars, double a_radius)
    {
        //don't do checks if this is the first step
        if (m_initial_var_arrays.size() > 0)
        {
            //first remind us which variable we're referring to
            auto alpha_new = a_vars[0];
            auto beta_new = a_vars[1];
            auto psi_new = a_vars[2];
            auto Psi_new = a_vars[3];
            auto psi_old = m_initial_var_arrays.back()[2];
            auto psi_central = m_initial_var_arrays[0][2];

            //if our solution has grown too large throw an exception
            if (std::abs(psi_new) > 2.0 * psi_central
            || std::abs(alpha_new) > 1.0e4 || std::abs(beta_new) > 1.0e4 ||
            std::abs(Psi_new) > 1.0e4)
            {
                throw std::runtime_error("Solution blow up. "\
                "Integration stopped.");
            }

            //if psi has changed sign, increase our count of the number of roots
            if (psi_new * psi_old < 0.0)
            {
                ++m_num_psi_roots;
            }

            //We only care about the ground state so stop integrating if the
            //number of roots exceeds 1
            if (m_num_psi_roots > 1)
            {
                throw std::runtime_error("This solution has more than 1 root"\
                " and is therefore not the ground state.");
            }
        }
        //finally store the result of this step and carry on
        m_initial_var_arrays.push_back(a_vars);
        m_initial_grid.push_back(a_radius);
    }
};

#endif /* BOSONSTARSOLUTIONOBSERVER_HPP_ */
