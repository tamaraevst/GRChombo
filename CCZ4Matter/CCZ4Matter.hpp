// Last edited K Clough 31.01.17

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

//!  Calculates RHS using CCZ4 including matter terms, and matter variable evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits from
     the CCZ4 class, which it uses to do the non matter evolution of variables.
     It then adds in the additional matter terms to the CCZ4 evolution (those including
     the stress energy tensor), and calculates the evolution of the matter variables.
     It does not assume a specific form of matter but is templated over a matter class
     matter_t. Please see the class SFMatter as an example of a matter_t.
     \sa CCZ4(), SFMatter()
*/

template <class matter_t>
class CCZ4Matter : public CCZ4 {

  //Use the variable definition in matter_t
  template<class data_t>
  using Vars=typename matter_t::template Vars<data_t>;
public:
  //!  Constructor of class CCZ4Matter
  /*!
       Inputs are the box driver and grid spacing, plus the CCZ4 evolution parameters and
       the matter parameters.
       It also takes the dissipation parameter sigma, and allows the formulation to be
       toggled between CCZ4 and BSSN. The default is CCZ4. It allows the user to set
       the value of Newton's constant, which is set to one by default.
  */
  CCZ4Matter(const FABDriverBase& driver, params_t params, const typename matter_t::matter_params_t matter_params,
             double dx, double sigma, int formulation = CCZ4::USE_CCZ4,
             double G_Newton = 1.0);

  //!  The compute member which calculates the RHS at each point in the box
  /*!
       \param ix the integer x coordinate of the current grid-cell.
       \param iy the integer y coordinate of the current grid-cell.
       \param iz the integer z coordinate of the current grid-cell.
       \return is void - the RHS is written directly to the grid in the function.
       \sa matter_rhs_equation()
  */
  template <class data_t>
  void compute(int ix, int iy, int iz);

protected:
  //! The function which calculates the RHS, given the vars and derivatives
  /*!
       \param vars the value of the variables at the point.
       \param d1 the value of the first derivatives of the variables.
       \param d2 the value of the second derivatives of the variables.
       \param advec the value of the advection terms beta^i d_i(var).
       \return is the RHS data for each variable at that point.
       \sa compute()
  */
  template <class data_t>
  void add_EMTensor_rhs(
      Vars<data_t> &matter_rhs,
      const Vars<data_t> &vars,
      const Vars< tensor<1,data_t> > &d1,
      const Vars< tensor<2,data_t> > &d2,
      const Vars<data_t> &advec);

  //! The matter params
  const typename matter_t::matter_params_t m_matter_params;

  //! Newton's constant, set to one by default.
  const double m_G_Newton;
};

#include "CCZ4Matter.impl.hpp"

#endif /* CCZ4MATTER_HPP_ */
