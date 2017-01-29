// Last edited K Clough 16.10.16

#if !defined(CONSTRAINTSMATTER_HPP_)
#error "This file should only be included through ConstraintsMatter.hpp"
#endif

#ifndef CONSTRAINTSMATTER_IMPL_HPP_
#define CONSTRAINTSMATTER_IMPL_HPP_

//TODO KClough: inline?
template <class matter_t>
ConstraintsMatter<matter_t>::ConstraintsMatter(
    const FABDriverBase& driver,
    const typename matter_t::matter_params_t matter_params,
    double dx,
    double G_Newton)
    : Constraints(driver, dx, 0.0), //Cosmological constant set to zero
      m_matter_params (matter_params),
      m_G_Newton (G_Newton) {}


template <class matter_t>
template <class data_t>
void ConstraintsMatter<matter_t>::compute(int x, int y, int z)
{
  //Calculate non matter contributions to Constraints
  vars_t<data_t> CCZ4_vars;
  m_driver.local_vars(CCZ4_vars);

  //Calculate first derivatives
  vars_t< tensor<1, data_t> > CCZ4_d1;
  FOR1(idir) m_deriv.diff1(CCZ4_d1, idir);

  //Calculate second derivatives
  vars_t< tensor<2,data_t> > CCZ4_d2;
  // Repeated derivatives
  FOR1(idir) m_deriv.diff2(CCZ4_d2, idir);
  // Mixed derivatives
  // Note: no need to symmetrise explicitely, this is done in mixed_diff2
  m_deriv.mixed_diff2(CCZ4_d2, 1, 0);
  m_deriv.mixed_diff2(CCZ4_d2, 2, 0);
  m_deriv.mixed_diff2(CCZ4_d2, 2, 1);

  // Get the non matter terms
  constraints_t<data_t> out = constraint_equations(CCZ4_vars, CCZ4_d1, CCZ4_d2);

  //Calculate EM Tensor and add matter terms, need advection and geometric objects
  //TODO K Clough: Once we template the Constraints class we won't need to calculate
  // d1 twice (or advec) as we can use the same vars_t object for both.
  matter_t my_matter(m_matter_params);
  typename matter_t::vars_t<data_t> vars;
  m_driver.local_vars(vars);

  typename matter_t::vars_t< tensor<1, data_t> > d1;
  FOR1(idir) m_deriv.diff1(d1, idir);

  typename matter_t::vars_t<data_t> advec;
  advec.assign(0.);
  FOR1(idir) m_deriv.add_advection(advec, vars.shift[idir], idir);

  // Inverse metric and Christoffel symbol
  auto h_UU = TensorAlgebra::compute_inverse(vars.h);
  auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

  // Energy Momentum tensor
  typename matter_t::emtensor_t<data_t> emtensor =
      my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL, advec);

  //Hamiltonain constraint
  out.Ham += -16.0*M_PI*m_G_Newton*emtensor.rho;

  //Momentum constraints
  FOR1(i)
  {
    out.Mom[i] += -8.0*M_PI*m_G_Newton*emtensor.Si[i];
  }

  //Write the rhs into the output FArrayBox
  m_driver.store_vars(out.Ham, c_Ham);
  m_driver.store_vars(out.Mom[0], c_Mom1);
  m_driver.store_vars(out.Mom[1], c_Mom2);
  m_driver.store_vars(out.Mom[2], c_Mom3);

}

#endif /* CONSTRAINTSMATTER_IMPL_HPP_ */
