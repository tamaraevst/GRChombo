/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_IMPL_HPP_
#define BOSONSTAR_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarBinarySearch.hpp" //for BosonStarBinarySearch class

inline BosonStar::BosonStar(BosonStar_params_t a_params_BosonStar,
                    Potential::params_t a_params_potential, double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential),
    m_1d_sol(a_params_BosonStar, a_params_potential)
{
}

void BosonStar::compute_1d_solution(const double a_max_radius)
{
    try
    {
        BosonStarBinarySearch<initial_data_t, initial_state_t> binary_search(
        m_params_BosonStar, m_params_potential);

        binary_search.shoot();
        auto sol = binary_search.getShootedSolution();

        m_1d_sol.makeFromPolarArealSolution(sol, a_max_radius);
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BosonStar::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    //conformal factor and lapse set to unity for now
    vars.chi = 1.;
    vars.lapse = 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
