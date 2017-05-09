// Last edited K Clough 08.05.17

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
#include "CCZ4Matter.hpp"
#include "CCZ4.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>

//!  Calculates the matter type specific elements such as the EMTensor and matter evolution
/*!
     This class is an example of a matter_t object which calculates the matter type specific
     elements for the RHS update and the evaluation of the constraints. This includes the
     Energy Momentum Tensor, and the matter evolution terms. In this case, a scalar field,
     the matter elements are phi and (minus) its conjugate momentum, Pi. It is templated over a
     potential function my_potential_t which the user must specify in a class. It assumes minimal 
     coupling of the field to gravity.
     \sa CCZ4Matter(), ConstraintsMatter()
*/
template <class my_potential_t>
class ScalarField
{
protected:
    //!The local copy of the potential
    my_potential_t my_potential;

public:
    //!  Constructor of class ScalarField, inputs are the matter parameters.
    ScalarField(const my_potential_t a_potential) : my_potential (a_potential) {}

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

    //! A structure for the potential data - the value of V and its gradient
    template <class data_t>
    struct potential_t {
        data_t  V_of_phi; //!< V(\phi)
        data_t  dVdphi; //!< Gradient of V(\phi)
    };

    //! The function which calculates the EM Tensor, given the vars and derivatives
    template <class data_t, template<typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,//!< the value of the variables at the point.
        const vars_t< tensor<1,data_t> >& d1,//!< the value of the first derivatives of the variables.
        const tensor<2, data_t>& h_UU,//!< the inverse metric (raised indices)
        const tensor<3, data_t>& chris_ULL,//!< the conformal chrisoffel symbol in ULL form.
        const vars_t<data_t> &advec//!< the value of the advection terms beta^i d_i(var)
    );

    //! The function which adds in the RHS for the matter field vars
    template <class data_t, template<typename> class vars_t>
    void add_matter_rhs(
        vars_t<data_t> &total_rhs, //!< contains the value of the RHS terms for all vars.
        const vars_t<data_t> &vars, //!< the value of the variables at the point.
        const vars_t< tensor<1,data_t> >& d1, //!< the value of the first derivatives of the variables.
        const vars_t< tensor<2,data_t> >& d2, //!< the value of the second derivatives of the variables.
        const vars_t<data_t> &advec); //!< the value of the advection terms beta^i d_i(var).

};

#include "ScalarField.impl.hpp"

#endif /* SCALARFIELD_HPP_ */
