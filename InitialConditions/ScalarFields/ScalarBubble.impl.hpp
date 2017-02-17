#if !defined(SCALARBUBBLE_HPP_)
#error "This file should only be included through ScalarBubble.hpp"
#endif

#ifndef SCALARBUBBLE_IMPL_HPP_
#define SCALARBUBBLE_IMPL_HPP_

#include "ScalarField.hpp"
#include "simd.hpp"

inline
ScalarBubble::ScalarBubble(const FABDriverBase& a_driver, ScalarField::params_t a_matter_params, double a_dx)
    : m_driver (a_driver), m_dx (a_dx), m_matter_params (a_matter_params)
{}

// Compute the value of the initial vars on the grid
template <class data_t>
void ScalarBubble::compute(int ix, int iy, int iz) {

  ScalarField::Vars<data_t> vars;
  vars.assign(0.); //Set only the non-zero components explicitly below
  Coordinates<data_t> coords(ix,iy,iz,m_dx);

  vars.phi = compute_phi(coords);
  vars.Pi = 0;

  vars.lapse = 1;

  vars.chi = 1;
  //Conformal metric is flat
  FOR1(i) vars.h[i][i] = 1.;

  m_driver.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t ScalarBubble::compute_phi(Coordinates<data_t> coords) {

  data_t rr2 =   pow(coords.x - m_matter_params.centerSF[0],2)
               + pow(coords.y - m_matter_params.centerSF[1],2)
               + pow(coords.z - m_matter_params.centerSF[2],2);

  double minimum_rr2 = 1e-12;
  auto r_is_too_small = simd_compare_lt(rr2, minimum_rr2);
  rr2 = simd_conditional(r_is_too_small, minimum_rr2, rr2);

  data_t rr = sqrt(rr2);
  double R0 = 5.0;
  data_t out_phi = m_matter_params.amplitudeSF*rr2*exp(-(rr-R0)*(rr-R0)/m_matter_params.widthSF);

  return out_phi;

}

#endif /* SCALARBUBBLE_IMPL_HPP_ */

