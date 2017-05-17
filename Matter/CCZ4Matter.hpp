// Last edited K Clough 16.02.17

#ifndef CCZ4MATTER_HPP_
#define CCZ4MATTER_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "VarsBase.hpp"
#include "CCZ4.hpp"
#include "Cell.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>

//!  Calculates RHS using CCZ4 including matter terms, and matter variable evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits from
     the CCZ4 class, which it uses to do the non matter evolution of variables.
     It then adds in the additional matter terms to the CCZ4 evolution (those including
     the stress energy tensor), and calculates the evolution of the matter variables.
     It does not assume a specific form of matter but is templated over a matter class
     matter_t. Please see the class ScalarField as an example of a matter_t.
     \sa CCZ4(), ScalarField()
*/

template <class matter_t>
class CCZ4Matter : public CCZ4
{
    //Use the variable definition in matter_t
    template<class data_t>
    using Vars=typename matter_t::template Vars<data_t>;

public:
    //!  Constructor of class CCZ4Matter
    /*!
       Inputs are the grid spacing, plus the CCZ4 evolution parameters and a matter object.
       It also takes the dissipation parameter sigma, and allows the formulation to be
       toggled between CCZ4 and BSSN. The default is CCZ4. It allows the user to set
       the value of Newton's constant, which is set to one by default.
    */
    CCZ4Matter(matter_t a_matter, params_t params,
               double dx, double sigma, int formulation = CCZ4::USE_CCZ4,
               double G_Newton = 1.0);

    //!  The compute member which calculates the RHS at each point in the box \sa matter_rhs_equation()
    template <class data_t>
    void compute(Cell current_cell);

protected:
    //! The function which adds in the EM Tensor terms to the CCZ4 rhs \sa compute()
    template <class data_t>
    void add_EMTensor_rhs(
        Vars<data_t> &matter_rhs, //!<the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!<the value of the variables at the point.
        const Vars< tensor<1,data_t> > &d1, //!<the value of the first derivatives of the variables.
        const Vars< tensor<2,data_t> > &d2, //!<the value of the second derivatives of the variables.
        const Vars<data_t> &advec //!<the value of the advection terms beta^i d_i(var).
    );

    // Class members
    matter_t my_matter;//!< The matter object, e.g. a scalar field.
    const double m_G_Newton;//!<Newton's constant, set to one by default.
};

#include "CCZ4Matter.impl.hpp"

#endif /* CCZ4MATTER_HPP_ */
