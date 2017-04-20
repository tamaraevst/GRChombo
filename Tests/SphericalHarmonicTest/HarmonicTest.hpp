#ifndef HARMONICTEST_HPP_
#define HARMONICTEST_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "FABDriverBase.hpp"
#include <vector>
#include <array>
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "Cell.hpp"

class HarmonicTest {
 public:
  HarmonicTest(const FABDriverBase& a_driver, IntVect a_center_vector, double a_dx)
        : m_driver (a_driver), m_dx (a_dx), m_center_vector (a_center_vector) {}

  template <class data_t>
  void compute(Cell current_cell);

 protected:
  const FABDriverBase& m_driver;
  double m_dx;
  IntVect m_center_vector;

  template<class data_t>
  data_t compute_harmonic(Coordinates<data_t> coords);

};

#include "HarmonicTest.impl.hpp" 

#endif /* HARMONICTEST_HPP_ */
