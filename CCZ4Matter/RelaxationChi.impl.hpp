// Last edited K Clough 16.10.16

#if !defined(RELAXATIONCHI_HPP_)
#error "This file should only be included through RelaxationChi.hpp"
#endif

#ifndef RELAXATIONCHI_IMPL_HPP_
#define RELAXATIONCHI_IMPL_HPP_

//TODO KClough: inline??
template <class matter_t>
RelaxationChi<matter_t>::RelaxationChi(
    const FABDriverBase& driver,
    SFMatter::matter_params_t matter_params,
    double dx,
    double relaxspeed,
    double G_Newton)
    : m_matter_params (matter_params),
      m_relaxspeed (relaxspeed),
      m_G_Newton (G_Newton),
      m_driver (driver),
      m_deriv (dx, m_driver) {}

template <class matter_t>
template <class data_t>
void RelaxationChi<matter_t>::compute(int ix, int iy, int iz) {

  //copy data from chombo gridpoint into local variables
  typename matter_t::vars_t<data_t> vars;
  m_driver.local_vars(vars);

  //work out first derivatives of variables on grid
  typename matter_t::vars_t< tensor<1, data_t> > d1;
  FOR1(idir) m_deriv.diff1(d1, idir);

  //work out second derivatives of variables on grid
  typename matter_t::vars_t< tensor<2,data_t> > d2;
  // Repeated derivatives
  FOR1(idir) m_deriv.diff2(d2, idir);
  // Mixed derivatives
  // Note: no need to symmetrise explicitely, this is done in mixed_diff2
  m_deriv.mixed_diff2(d2, 1, 0);
  m_deriv.mixed_diff2(d2, 2, 0);
  m_deriv.mixed_diff2(d2, 2, 1);

  // Calculate advection components
  typename matter_t::vars_t<data_t> advec;
  advec.assign(0.);
  FOR1(idir) m_deriv.add_advection(advec, vars.shift[idir], idir);

  //work out RHS including advection
  typename matter_t::vars_t<data_t> rhs = rhs_equation(vars, d1, d2, advec);

  //    No dissipation in relaxation for now but may add it
  //    FOR1(idir) m_deriv.add_dissipation(rhs, m_sigma, idir);

  //Write the rhs into the output FArrayBox
  m_driver.store_vars(rhs);
}

template <class matter_t>
template <class data_t>
typename matter_t::vars_t<data_t> RelaxationChi<matter_t>::rhs_equation(
    const typename matter_t::vars_t<data_t> &vars,
    const typename matter_t::vars_t< tensor<1,data_t> >& d1,
    const typename matter_t::vars_t< tensor<2,data_t> >& d2,
    const typename matter_t::vars_t<data_t> &advec) {

  typename matter_t::vars_t<data_t> rhs;
  rhs.assign(0);

  using namespace TensorAlgebra;

  auto h_UU = compute_inverse(vars.h);
  auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

  //Calculate elements of the decomposed stress energy tensor and ricci tensor
  matter_t my_matter(m_matter_params);
  auto emtensor =  my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL, advec);
  auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);
  auto A_UU       = TensorAlgebra::raise_all(vars.A, h_UU);
  data_t tr_AA    = TensorAlgebra::compute_trace(vars.A, A_UU);

  //Calculate the relaxation RHS for chi, all other vars RHS zero
  //Could have called ConstraintsMatter here, but it is hardly worth it
  rhs.chi =  m_relaxspeed*(ricci.scalar+(GR_SPACEDIM-1.)*vars.K*vars.K/GR_SPACEDIM
                           - tr_AA - 16.0*M_PI*m_G_Newton*emtensor.rho);

  return rhs;
}

#endif /* RELAXATIONCHI_IMPL_HPP_ */
