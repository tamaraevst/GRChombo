#ifndef CCZ4MATTER_HPP_
#define CCZ4MATTER_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "SFMatter.hpp"
#include "VarsBase.hpp"
#include "CCZ4.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>

template <class matter_t>
class CCZ4Matter : public CCZ4 {
 public:
  CCZ4Matter(const FABDriverBase& driver, params_t params, double dx,
             double sigma, int formulation = CCZ4::USE_CCZ4,
             double G_Newton = 1.0);

  template <class data_t>
  void compute(int ix, int iy, int iz);

 protected:
  template <class data_t>
  typename matter_t::vars_t<data_t> matter_rhs_equation(
      const typename matter_t::vars_t<data_t> &vars,
      const typename matter_t::vars_t< tensor<1,data_t> > &d1,
      const typename matter_t::vars_t< tensor<2,data_t> > &d2,
      const typename matter_t::vars_t<data_t> &advec);

  double m_G_Newton;
};

#include "CCZ4Matter.impl.hpp"

#endif /* CCZ4MATTER_HPP_ */
