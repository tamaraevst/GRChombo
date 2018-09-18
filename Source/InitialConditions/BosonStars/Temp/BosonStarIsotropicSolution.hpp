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

//! This class constructs grid values for the boson star solution isotropic
//! coordinates given a rescaled solution in polar-areal coordinates

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarIsotropicSolution
{
public:
    //! Constructor
    BosonStarIsotropicSolution(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution, BosonStar::params_t a_params_BosonStar,
        Potential::params_t a_params_potential, const double a_max_radius);

    //! This function calculates the isotropic grid from the polar areal grid.
    void calculateIsotropicGrid(const double a_max_radius);

private:
    BosonStar::params_t m_params_BosonStar;
    Potential::params_t m_params_potential;
    BosonStarSolution<initial_data_t, initial_state_t>& m_polar_areal_solution;
    initial_data_t<double> m_isotropic_grid = {};
    initial_data_t<double> m_polar_areal_grid = {};
    initial_data_t<double> m_chi_array = {};
    initial_data_t<double> m_lapse_array = {};
    initial_data_t<double> m_phi_array = {};
};

#include "BosonStarIsotropicSolution.impl.hpp"

#endif /* BOSONSTARISOTROPICSOLUTION_HPP_ */
