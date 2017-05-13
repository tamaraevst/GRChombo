#if !defined(SCALARBUBBLE_HPP_)
#error "This file should only be included through ScalarBubble.hpp"
#endif

#ifndef SCALARBUBBLE_IMPL_HPP_
#define SCALARBUBBLE_IMPL_HPP_

inline
ScalarBubble::ScalarBubble(const FABDriverBase& a_driver, params_t a_params, double a_dx)
    : m_driver (a_driver), m_dx (a_dx), m_params (a_params)
{}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarBubble::compute(Cell current_cell)
{
    ScalarField<Potential>::Vars<data_t> vars;
    vars.assign(0.); //Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell,m_dx);

    //set the field vars
    vars.phi = compute_phi(coords);
    vars.Pi = 0;

    //start with unit lapse and flat metric (must be relaxed for chi)
    vars.lapse = 1;
    vars.chi = 1;

    //conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    //Store the initial values of the variables
    m_driver.store_vars(vars, current_cell);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarBubble::compute_phi(Coordinates<data_t> coords)
{
    data_t rr = coords.get_radius(m_params.centerSF);
    data_t rr2 = rr*rr;
    data_t out_phi = m_params.amplitudeSF*rr2*exp(-pow(rr-m_params.r_zero/m_params.widthSF, 2.0));

    return out_phi;
}

#endif /* SCALARBUBBLE_IMPL_HPP_ */

