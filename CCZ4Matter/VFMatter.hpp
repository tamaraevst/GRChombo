// Last edited K Clough 30.01.2017

#ifndef VFMATTER_HPP_
#define VFMATTER_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "VarsBase.hpp"
#include "CCZ4.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>

//!  Calculates the matter type specific elements such as the EMTensor and matter evolution
/*!
     This class is an example of a matter_t object which calculates the matter type specific
     elements for the RHS update and the evaluation of the constraints. This includes the
     Energy Momentum Tensor, and the matter evolution terms. In this case, a vector field,
     the matter elements are phi, A_i and E_i. It assumes minimal coupling of the field to gravity.
     \sa CCZ4Matter(), ConstraintsMatter()
*/

class VFMatter {
 public:
  //! A structure for the input params for VF properties and initial conditions
  struct matter_params_t {
      double amplitudeSF;
      std::vector<double> centerSF;
      double widthSF;
      double vector_mass;
  };

 protected:
  //! The local copy of the matter params
  matter_params_t m_matter_params;

 public:
  //!  Constructor of class VFMatter
  /*!
       Inputs are the matter parameters.
  */
  VFMatter(matter_params_t matter_params) : m_matter_params (matter_params) {}

  // May not be needed after templating vars_t
  template <class data_t>
  struct vars_t : VarsBase<data_t> {

    using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later
    using VarsBase<data_t>::define_symmetric_enum_mapping;

    data_t chi;
    tensor<2, data_t> h;
    data_t K;
    tensor<2, data_t> A;
    tensor<1, data_t> Gamma;
    data_t Theta;
    data_t lapse;
    tensor<1, data_t> shift;
    tensor<1, data_t> B;
    data_t phi;

    vars_t();
  };

  //! A structure for the decomposed elements of the Energy Momentum Tensor in 3+1D
  template <class data_t>
  struct emtensor_t {
      tensor<2, data_t> Sij; // S_ij = T_ij
      tensor<1,data_t>  Si; // S_i = T_ia_n^a
      data_t            S; // S = S^i_i
      data_t            rho; // rho = T_ab n^a n^b
  };

  //! A structure for the potential data - the value of V and its gradient
  template <class data_t>
  struct potential_t {
      data_t  V_of_phi; //V(\phi)
      data_t  dVdphi; // Gradient of V(\phi)
  };

  //! The function which calculates the EM Tensor, given the vars and derivatives
  /*!
      \param vars the value of the variables at the point.
      \param d1 the value of the first derivatives of the variables.
      \param h_UU the inverse metric (raised indices)
      \param chris_ULL is the conformal chrisoffel symbol in ULL form
      \param advec the value of the advection terms beta^i d_i(var).
      \return is the EM Tensor structure
      \sa compute_potential()
  */
  template <class data_t>
  emtensor_t<data_t> compute_emtensor(
      const vars_t<data_t> &vars,
      const vars_t< tensor<1,data_t> >& d1,
      const tensor<2, data_t>& h_UU,
      const tensor<3, data_t>& chris_ULL,
      const vars_t<data_t> &advec
  );

  //! The function which calculates the potential function, given the field value
  /*!
      \param phi is the value of the field at the point
      \return is the potential structure - the value of V and its derivative
      \sa compute_emtensor()
  */
  template <class data_t>
  potential_t<data_t> compute_potential(const data_t phi);

  //! The function which calculates the total RHS, given the CCZ4 RHSs, vars and derivatives
  /*!
      \param CCZ4_rhs the value of the RHS terms calculated by the CCZ4 class.
      \param vars matter_rhs the value of the RHS terms in the CCZ4 evolution due to matter.
      \param vars the value of the variables at the point.
      \param d1 the value of the first derivatives of the variables.
      \param d2 the value of the second derivatives of the variables.
      \param advec the value of the advection terms beta^i d_i(var).
      \return is the RHS data for each variable at that point.
      \sa compute_potential()
  */
  template <class data_t>
  vars_t<data_t> compute_total_rhs(
      const CCZ4::vars_t<data_t> &CCZ4_rhs,
      const vars_t<data_t> &matter_rhs,
      const vars_t<data_t> &vars,
      const vars_t< tensor<1,data_t> >& d1,
      const vars_t< tensor<2,data_t> >& d2,
      const vars_t<data_t> &advec);

};

#include "VFMatter.impl.hpp"

#endif /* VFMATTER_HPP_ */
