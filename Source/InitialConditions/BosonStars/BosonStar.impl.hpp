/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_IMPL_HPP_
#define BOSONSTAR_IMPL_HPP_

inline BosonStar::BosonStar(params_t a_params, double a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BosonStar::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    //Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

    // conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */
