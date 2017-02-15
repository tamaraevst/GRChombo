#ifndef BUBBLESF_HPP_
#define BUBBLESF_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "FABDriverBase.hpp"
#include <vector>
#include "tensor.hpp"
#include <array>
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "SFMatter.hpp"

class BubbleSF {
 public:
  BubbleSF(const FABDriverBase& a_driver, SFMatter::matter_params_t a_matter_params, double a_dx);

  template <class data_t>
  void compute(int ix, int iy, int iz);

 protected:
  const FABDriverBase& m_driver;
  double m_dx;

  template<class data_t>
  data_t compute_phi(Coordinates<data_t> coords);

 public:
  const SFMatter::matter_params_t m_matter_params;
};

#include "BubbleSF.impl.hpp"

#endif /* BUBBLESF_HPP_ */
