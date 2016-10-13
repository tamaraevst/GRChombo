#ifndef RELAXATIONCHI_HPP_
#define RELAXATIONCHI_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "SFMatter.hpp"
#include "VarsBase.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>

template <class matter_t>
class RelaxationChi {
 public:
  RelaxationChi(const FABDriverBase& driver, double dx, double relaxspeed,
    double G_Newton = 1.0);

    template <class data_t>
    void compute(int ix, int iy, int iz);

 protected:
  const double m_relaxspeed;
  const double m_G_Newton;
  const FABDriverBase& m_driver;
  const FourthOrderDerivatives m_deriv;

  template <class data_t>
  typename matter_t::vars_t<data_t> rhs_equation(
      const typename matter_t::vars_t<data_t> &vars,
      const typename matter_t::vars_t< tensor<1,data_t> > &d1,
      const typename matter_t::vars_t< tensor<2,data_t> > &d2,
      const typename matter_t::vars_t<data_t> &advec);

};

#include "RelaxationChi.impl.hpp"

#endif /* RELAXATIONCHI_HPP_ */
