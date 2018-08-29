/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARINTEGRATOR_HPP_
#define BOSONSTARINTEGRATOR_HPP_

#include "BosonStar.hpp"
#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "BosonStarSolution.hpp"
#include <boost/numeric/odeint.hpp>


//! Class that uses odeint to integrate the equations for a static, spherically
//! symmetric boson star.
template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarIntegrator
{
    typedef initial_data_t<double> initial_grid_t;

public:
    //! Constructor
    BosonStarIntegrator(BosonStar::params_t a_params_BosonStar,
        Potential::params_t a_params_potential);

    //! Clear internally arrays. Called by doIntegration()
    void clearArrays();

    //! Do integration
    void doIntegration(const double a_alpha_central);

    //! Returns the solution in a BosonStarSolution object.
    BosonStarSolution<initial_data_t, initial_state_t> getSolution();

protected:
    BosonStar::params_t m_params_BosonStar;
    BosonStarRHS m_boson_star_rhs;
    initial_data_t<initial_state_t> m_initial_var_arrays;
    initial_grid_t m_initial_grid;
    BosonStarSolutionObserver<initial_data_t, initial_state_t> m_sol_observer;
};

#include "BosonStarIntegrator.impl.hpp"

#endif /* BOSONSTARINTEGRATOR_HPP_ */
