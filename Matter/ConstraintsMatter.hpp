// Last edited K Clough 08.05.17

#ifndef CONSTRAINTSMATTER_HPP_
#define CONSTRAINTSMATTER_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FourthOrderDerivatives.hpp"
#include "CCZ4Geometry.hpp"
#include "Constraints.hpp"
#include <array>
#include "Cell.hpp"

//!  Calculates the Hamiltonain and Momentum constraints with matter fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point in a
     box. It inherits from the Constraints class which calculates the constraints without
     the matter terms. It adds in the matter terms for a given matter class matter_t, which
     must provide it with the Energy Momentum tensor. For an example of a matter_t class
     see ScalarField.
     \sa Constraints(), ScalarField()
*/
template <class matter_t>
class ConstraintsMatter : public Constraints
{
    //Use the variable definition in matter_t
    template<class data_t>
    using Vars=typename matter_t::template Vars<data_t>;

public:
    //!  Constructor of class ConstraintsMatter
    /*!
         Takes in the grid spacing, and matter object plus
         optionally the value of Newton's constant, which is set to one by default.
    */
    ConstraintsMatter(const matter_t a_matter,
                      double dx, double G_Newton = 1.0);

    //! The compute member which calculates the constraints at each point in the box
    template <class data_t>
    void compute(Cell current_cell);

protected:
    matter_t my_matter; //!< The matter object, e.g. a scalar field
    double m_G_Newton; //!< Newton's constant, set to one by default.

};

#include "ConstraintsMatter.impl.hpp"

#endif /* CONSTRAINTSMATTER_HPP_ */
