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


//!  Calculates RHS for relaxation of the conformal factor, for initial conditions
/*!
     The class calculates the RHS evolution for the relaxation of the conformal
     factor chi, in order to satisfy the Hamiltonian constraint for some general
     matter configuration. It assumes that the other variables (K, h_ij etc) are
     fixed and satisfy the Momentum constraint. It calculates the RHS at each step
     as the relaxation speed multiplied by the error in the Hamiltonian constraint.
     It is extremely inefficient and takes a long time to converge, but it works.
     \sa m_relaxspeed()
*/

template <class matter_t>
class RelaxationChi {
 public:
  //! Constructor of class RelaxationChi
  /*!
      Takes in the box driver and the grid spacing, plus the relaxation speed and
      value of Newton's constant, which is set to one by default.
  */
  RelaxationChi(const FABDriverBase& driver, double dx, double relaxspeed,
    double G_Newton = 1.0);

  //! The compute member which calculates the RHS at each point in the box
  /*!
      \param ix the x index of the point on the level.
      \param iy the y index of the point on the level.
      \param iz the z index of the point on the level.
      \return is void - the RHS is written directly to the grid in the function
      \sa rhs_equation()
  */
  template <class data_t>
  void compute(int ix, int iy, int iz);

 protected:
  //! The coefficient of the Hamiltonian used to determine relaxation speed.
  const double m_relaxspeed;
  //! Newton's constant, set to one by default.
  const double m_G_Newton;
  //! The driver for the array box
  const FABDriverBase& m_driver;
  //! An object for calculating derivatives of the variables at the point
  const FourthOrderDerivatives m_deriv;

  //! The function which calculates the RHS, given the vars and derivatives
  /*!
      \param vars the value of the variables at the point.
      \param d1 the value of the first derivatives of the variables.
      \param d2 the value of the second derivatives of the variables.
      \param advec the value of the advection terms beta^i d_i(var)
      \return is the RHS data for each variable at that point
      \sa compute()
  */
  template <class data_t>
  typename matter_t::vars_t<data_t> rhs_equation(
      const typename matter_t::vars_t<data_t> &vars,
      const typename matter_t::vars_t< tensor<1,data_t> > &d1,
      const typename matter_t::vars_t< tensor<2,data_t> > &d2,
      const typename matter_t::vars_t<data_t> &advec);

};

#include "RelaxationChi.impl.hpp"

#endif /* RELAXATIONCHI_HPP_ */
