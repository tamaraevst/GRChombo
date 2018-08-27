/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARBINARYSEARCH_HPP_
#define BOSONSTARBINARYSEARCH_HPP_

#include "BosonStar.hpp"
#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "ComplexPotential"
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <iostream>
#include <stdexcept>

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarBinarySearch
{
public:
    //! Constructor 1 which does not require computed solutions
    BosonStarBinarySearch(BosonStar::params_t a_params_CSF,
        Potential::params_t a_params_potential, double alpha_central_min,
        double alpha_central_max);

    //! Constructor 2 which requires pre-computed solutions
    BosonStarBinarySearch(BosonStar::params_t a_params_CSF,
        Potential::params_t a_params_potential, double alpha_central_min,
        double alpha_central_max, initial_data_t<initial_state_t> sol_min,
        initial_data_t<initial_state_t> sol_max);

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
    void shoot(const double a_tol);

    /*! This function gets the solution and returns it to the caller in the
    variable out.
    */
    void getSolution(initial_data_t<initial_state_t> &out) const;

    /*! This function allows the user to change the binary search interval for
    alpha_central
    */
    void setInterval(double alpha_central_min, double alpha_central_max);

private:
    BosonStar::params_t m_params_CSF;
    Potential::params_t m_params_potential;
    double m_alpha_central_min;
    double m_alpha_central_max;
    double m_alpha_central_mid;

    initial_data_t<initial_state_t> m_sol_min;
    initial_data_t<initial_state_t> m_sol_max;
    initial_data_t<initial_state_t> m_sol_mid;
};

#include "BosonStarBinarySearch.impl.hpp"

#endif /* BOSONSTARBINARYSEARCH_HPP_ */
