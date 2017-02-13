// Last edited K Clough 31.01.17

#ifndef SFMATTER_HPP_
#define SFMATTER_HPP_

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
     Energy Momentum Tensor, and the matter evolution terms. In this case, a scalar field,
     the matter elements are phi and (minus) its conjugate momentum, Pi. It also specifies
     the form of the potential function. It assumes minimal coupling of the field to gravity.
     \sa CCZ4Matter(), ConstraintsMatter()
*/

class SFMatter {
 public:
  //! A structure for the input params for SF properties and initial conditions
  struct matter_params_t {
      double amplitudeSF;
      std::vector<double> centerSF;
      double widthSF;
      double scalar_mass;
  };

 protected:
  //! The local copy of the matter params
  matter_params_t m_matter_params;

 public:
  //!  Constructor of class SFMatter
  /*!
       Inputs are the matter parameters.
  */
  SFMatter(matter_params_t matter_params) : m_matter_params (matter_params) {}

  //! Structure containing all the rhs variables for the gravity and matter fields
  template <class data_t>
  struct Vars : VarsBase<data_t> {

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
    data_t Pi;

    Vars();
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
      const Vars<data_t> &vars,
      const Vars< tensor<1,data_t> >& d1,
      const tensor<2, data_t>& h_UU,
      const tensor<3, data_t>& chris_ULL,
      const Vars<data_t> &advec
  );

  //! The function which calculates the potential function, given the field value
  /*!
      \param phi is the value of the field at the point
      \return is the potential structure - the value of V and its derivative
      \sa compute_emtensor()
  */
  template <class data_t>
  potential_t<data_t> compute_potential(const data_t phi);

  //! The function which adds in the matter field RHS, given the vars and derivatives
  /*!
      \param matter_rhs contains the value of the RHS terms for all vars.
      \param vars the value of the variables at the point.
      \param d1 the value of the first derivatives of the variables.
      \param d2 the value of the second derivatives of the variables.
      \param advec the value of the advection terms beta^i d_i(var).
      \sa compute_potential()
  */
  template <class data_t>
  void add_matter_rhs(
      Vars<data_t> &total_rhs,
      const Vars<data_t> &vars,
      const Vars< tensor<1,data_t> >& d1,
      const Vars< tensor<2,data_t> >& d2,
      const Vars<data_t> &advec);

};

#include "SFMatter.impl.hpp"

#endif /* SFMATTER_HPP_ */
