// Last edited K Clough 31.01.17

#ifndef SCALARFIELD_HPP_
#define SCALARFIELD_HPP_

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

class ScalarField {
 public:
  //! A structure for the input params for scalar field properties and initial conditions
  struct params_t {
      double amplitudeSF;//!< Amplitude of bump in initial SF bubble
      std::vector<double> centerSF;//!< Centre of bump in initial SF bubble
      double widthSF;//!< Width of bump in initial SF bubble
      double scalar_mass;//!< Mass of the scalar field
  };

 protected:
  params_t m_params;//!< The local copy of the matter params

 public:
  //!  Constructor of class ScalarField, inputs are the matter parameters.
  ScalarField(params_t params) : m_params (params) {}

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
      tensor<2, data_t> Sij; //!< S_ij = T_ij
      tensor<1,data_t>  Si; //!< S_i = T_ia_n^a
      data_t            S; //!< S = S^i_i
      data_t            rho; //!< rho = T_ab n^a n^b
  };

  //! A structure for the potential data - the value of V and its gradient
  template <class data_t>
  struct potential_t {
      data_t  V_of_phi; //!< V(\phi)
      data_t  dVdphi; //!< Gradient of V(\phi)
  };

  //! The function which calculates the EM Tensor, given the vars and derivatives \sa compute_potential()
  template <class data_t>
  emtensor_t<data_t> compute_emtensor(
      const Vars<data_t> &vars,//!< the value of the variables at the point.
      const Vars< tensor<1,data_t> >& d1,//!< the value of the first derivatives of the variables.
      const tensor<2, data_t>& h_UU,//!< the inverse metric (raised indices)
      const tensor<3, data_t>& chris_ULL,//!< the conformal chrisoffel symbol in ULL form.
      const Vars<data_t> &advec//!< the value of the advection terms beta^i d_i(var)
  );

  //! The function which calculates the potential function, given the field value \sa compute_emtensor()
  template <class data_t>
  potential_t<data_t> compute_potential(const data_t phi);

  //! The function which adds in the RHS for the matter field vars \sa compute_potential()
  template <class data_t>
  void add_matter_rhs(
      Vars<data_t> &total_rhs, //!< contains the value of the RHS terms for all vars.
      const Vars<data_t> &vars, //!< the value of the variables at the point.
      const Vars< tensor<1,data_t> >& d1, //!< the value of the first derivatives of the variables.
      const Vars< tensor<2,data_t> >& d2, //!< the value of the second derivatives of the variables.
      const Vars<data_t> &advec); //!< the value of the advection terms beta^i d_i(var).

};

#include "ScalarField.impl.hpp"

#endif /* SCALARFIELD_HPP_ */
