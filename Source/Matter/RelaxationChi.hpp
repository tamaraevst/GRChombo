// Last edited K Clough 17.05.17

#ifndef RELAXATIONCHI_HPP_
#define RELAXATIONCHI_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Matter.hpp"
#include "VarsTools.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include <array>
#include "Cell.hpp"

//!  Calculates RHS for relaxation of the conformal factor, for initial conditions
/*!
     The class calculates the RHS evolution for the relaxation of the conformal
     factor chi, in order to satisfy the Hamiltonian constraint for some general
     matter configuration. It assumes that the other variables (K, h_ij etc) are
     fixed and satisfy the Momentum constraint. It calculates the RHS at each step
     as the relaxation speed multiplied by the error in the Hamiltonian constraint.
     It is extremely inefficient and takes a long time to converge, but it works
     provided that there is indeed a sensible solution (note that it can converge
     to \chi = 0 everywhere if there is not, so one must always check it is converging
     on something sensible). Note that the relaxation speed should be set to
     a value less than 2/15*dx_min/courant_factor for numerical stability.
     \sa m_relax_speed()
*/

template <class matter_t>
class RelaxationChi {

    //Use the variable definitions in CCZ4Matter
    template <class data_t>
    using Vars = typename CCZ4Matter<matter_t>::template Vars<data_t>;

    template<class data_t>
    using Diff2Vars = typename CCZ4Matter<matter_t>::template Diff2Vars<data_t>;

public:
    //! Constructor of class RelaxationChi
    /*!
        Takes in the grid spacing, plus the relaxation speed, a matter object
        and the value of Newton's constant, which is set to one by default.
    */
    RelaxationChi(matter_t a_matter, double dx, double relax_speed, double G_Newton = 1.0);

    //! The compute member which calculates the RHS at each point in the box \sa rhs_equation()
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

protected:
    matter_t my_matter;//!< The matter object, e.g. a scalar field.
    const double m_relax_speed;//!< The coefficient of the Hamiltonian used to set relaxation speed.
    const double m_G_Newton;//!< Newton's constant, set to one by default.
    const FourthOrderDerivatives m_deriv;//!< An object for calculating derivatives of the variables

    //! The function which calculates the RHS, given the vars and derivatives \sa compute()
    template <class data_t>
    void rhs_equation(
        Vars<data_t> &rhs,//!<the RHS data for each variable at that point.
        const Vars<data_t> &vars,//!< the value of the variables at the point.
        const Vars< tensor<1,data_t> > &d1,//!< the value of the first derivatives of the variables.
        const Diff2Vars< tensor<2,data_t> > &d2,//!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec //!< advec the value of the advection terms beta^i d_i(var)
        )
    const;

};

#include "RelaxationChi.impl.hpp"

#endif /* RELAXATIONCHI_HPP_ */
