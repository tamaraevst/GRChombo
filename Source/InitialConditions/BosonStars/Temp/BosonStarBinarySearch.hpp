/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARBINARYSEARCH_HPP_
#define BOSONSTARBINARYSEARCH_HPP_

#include "BosonStar.hpp"
#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "BosonStarSolution.hpp"
#include "ComplexPotential.hpp"
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarBinarySearch
{
    typedef initial_data_t<double> initial_grid_t;

public:
    //! Constructor
    BosonStarBinarySearch(BosonStar::params_t a_params_CSF,
        Potential::params_t a_params_potential,
        BosonStarSolution<initial_data_t, initial_state_t> &a_sol_min,
        BosonStarSolution<initial_data_t, initial_state_t> &a_sol_max);

    /*! Function called by the constructor to check the inputted values of
    alpha_central_min and alpha_central_max bound the desired value of
    alpha_central. This checks that the sign of psi at the edge of the
    computational domain differs between the two solutions and that they have
    the correct number of roots (0 for solution corresponding to
    alpha_central_max and 1 for solution corresponding to alpha_central min) so
    that we are getting the ground state
    */
    bool checkValidInterval();

    /*! This function performs the binary search shooting algorithm to find a
    value of alpha that is within a_tol of the true value
    */
    void shoot(const double a_tol, const double a_max_radius = 100.0);

    /*! This function gets the solution and returns it to the caller in the
    variable out.
    */
    BosonStarSolution<initial_data_t, initial_state_t>& getSolution();

private:
    BosonStar::params_t m_params_CSF;
    Potential::params_t m_params_potential;
    double m_alpha_central_min;
    double m_alpha_central_max;
    double m_alpha_central_mid;

    BosonStarSolution<initial_data_t, initial_state_t> m_sol_min;
    BosonStarSolution<initial_data_t, initial_state_t> m_sol_max;
    BosonStarSolution<initial_data_t, initial_state_t> m_sol_mid;
};

#include "BosonStarBinarySearch.impl.hpp"

#endif /* BOSONSTARBINARYSEARCH_HPP_ */
