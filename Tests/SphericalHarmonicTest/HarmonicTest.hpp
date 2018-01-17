#ifndef HARMONICTEST_HPP_
#define HARMONICTEST_HPP_

#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include "tensor.hpp"

class HarmonicTest
{
  public:
    HarmonicTest(std::array<double, CH_SPACEDIM> a_center_vector, double a_dx)
        : m_dx(a_dx), m_center_vector(a_center_vector)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    std::array<double, CH_SPACEDIM> m_center_vector;

    template <class data_t>
    data_t compute_harmonic(Coordinates<data_t> coords) const;
};

#include "HarmonicTest.impl.hpp"

#endif /* HARMONICTEST_HPP_ */
