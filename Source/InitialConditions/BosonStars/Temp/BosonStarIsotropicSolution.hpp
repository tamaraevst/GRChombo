/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARISOTROPICSOLUTION_HPP_
#define BOSONSTARISOTROPICSOLUTION_HPP_

#include "BosonStarSolution.hpp"
#include "ComplexPotential.hpp"
#include <cmath>
#include <algorithm>
#include "SplineInterpolator.hpp"
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <limits>

//! This class constructs interpolation functions for the boson star solution
//! in isotropic coordinates given a rescaled solution in polar-areal coordinates

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarIsotropicSolution
{
public:
    tools::spline<initial_data_t> m_chi;
    tools::spline<initial_data_t> m_lapse;
    tools::spline<initial_data_t> m_phi;

    //! Constructor
    BosonStarIsotropicSolution(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution, BosonStar::params_t a_params_BosonStar,
        Potential::params_t a_params_potential, const double a_max_radius,
        const double a_G_Newton = 1.0);

    //! This function calculates the isotropic grid from the polar areal grid.
    void calculateIsotropicGrid(const double a_max_radius);

    //! Adds the points for the conformal factor chi interpolation function
    void construct_chi();

    //! Adds the points for the lapse and scalar field interpolation functions
    void construct_phi_and_lapse();


private:
    BosonStar::params_t m_params_BosonStar;
    Potential::params_t m_params_potential;
    BosonStarSolution<initial_data_t, initial_state_t>& m_polar_areal_solution;
    initial_data_t<double> m_isotropic_grid = {};
    initial_data_t<double> m_polar_areal_grid = {};
    const double m_G_Newton;
};

#include "BosonStarIsotropicSolution.impl.hpp"

#endif /* BOSONSTARISOTROPICSOLUTION_HPP_ */
