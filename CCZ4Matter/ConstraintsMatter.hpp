//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTSMATTER_HPP_
#define CONSTRAINTSMATTER_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "CCZ4Geometry.hpp"
#include "SFMatter.hpp"
#include "Constraints.hpp"
#include <array>

template <class matter_t>
class ConstraintsMatter : public Constraints {
 public:
  ConstraintsMatter(const FABDriverBase& driver, double dx, double G_Newton);

  template <class data_t>
  void compute(int x, int y, int z);

 protected:
  double m_G_Newton;

};

#include "ConstraintsMatter.impl.hpp"

#endif /* CONSTRAINTSMATTER_HPP_ */
