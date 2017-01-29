#if !defined(BUBBLESF_HPP_)
#error "This file should only be included through BubbleSF.hpp
#endif

#ifndef BUBBLESF_IMPL_HPP_
#define BUBBLESF_IMPL_HPP_

#include "BubbleSF.hpp"
#include "SFMatter.hpp"
#include "simd.hpp"

template <class data_t>
void BubbleSF::compute(int ix, int iy, int iz) {

  SFMatter::vars_t<data_t> vars;
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

template <class data_t>
data_t BubbleSF::compute_phi(Coordinates<data_t> coords) {

  data_t rr2 = pow(coords.x - m_matter_params.centerSF[0],(decltype(coords.x))2)
                 + pow(coords.y - m_matter_params.centerSF[1],2)
                 + pow(coords.z - m_matter_params.centerSF[2],2);

  double minimum_rr2 = 1e-12;
  auto r_is_too_small = simd_compare_lt(rr2, minimum_rr2);
  rr2 = simd_conditional(r_is_too_small, minimum_rr2, rr2);

  data_t rr = sqrt(rr2);
  double R0 = 5.0;
  data_t phiout = m_matter_params.amplitudeSF*rr2*exp(-(rr-R0)*(rr-R0)/m_matter_params.widthSF);

  return phiout;

}

#endif /* BUBBLESF_IMPL_HPP_ */

