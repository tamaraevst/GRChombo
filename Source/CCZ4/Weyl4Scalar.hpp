/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYL4SCALAR_HPP_
#define WEYL4SCALAR_HPP_

#include "BSSNVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include <array>

//! Struct for the E and B fields
template <class data_t> struct EBFields_t
{
    Tensor<2, data_t> E; //!< Electric component of Weyltensor
    Tensor<2, data_t> B; //!< Magnetic component of Weyltensor
};

//! Struct for the null tetrad
template <class data_t> struct Tetrad_t
{
    Tensor<1, data_t> u; //!< the vector u^i
    Tensor<1, data_t> v; //!< the vector v^i
    Tensor<1, data_t> w; //!< the vector w^i
};

//! Struct for the Newman Penrose scalar
template <class data_t> struct NPScalar_t
{
    data_t Real; // Real component
    data_t Im;   // Imaginary component
};

//!  Calculates the Weyl4 scalar for spacetimes without matter content
/*!
   This class calculates the Weyl4 scalar real and im parts using definitions
   from Alcubierres book "Introduction to 3+1 Numerical Relativity". We use a
   decomposition of the Weyl tensor in electric and magnetic parts, which then
   is used together with the tetrads defined in "gr-qc/0104063" to calculate the
   Weyl4 scalar.
*/
class Weyl4Scalar
{
  public:
    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = BSSNVars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;

    //! Constructor of class Weyl4Scalar
    /*!
        Takes in the centre for the calculation of the tetrads, and grid spacing
    */
    Weyl4Scalar(const std::array<double, CH_SPACEDIM> a_center,
                const double a_dx)
        : m_center(a_center), m_dx(a_dx), m_deriv(a_dx)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars
    const double m_dx;                    //!< the grid spacing

    //! Calculation of Weyl_4 scalar
    template <class data_t>
    NPScalar_t<data_t> compute_Weyl4(const EBFields_t<data_t> &ebfields,
                                     const Vars<data_t> &vars,
                                     const Vars<Tensor<1, data_t>> &d1,
                                     const Diff2Vars<Tensor<2, data_t>> &d2,
                                     const Coordinates<data_t> &coords) const;

    //! Calculation of the tetrads
    template <class data_t>
    Tetrad_t<data_t>
    compute_null_tetrad(const Vars<data_t> &vars,
                        const Coordinates<data_t> &coords) const;

    //! Calulation of the decomposition of the Weyl tensor in Electric and
    //! Magnetic fields
    template <class data_t>
    EBFields_t<data_t>
    compute_EB_fields(const Vars<data_t> &vars,
                      const Vars<Tensor<1, data_t>> &d1,
                      const Diff2Vars<Tensor<2, data_t>> &d2,
                      const Coordinates<data_t> &coords) const;
};

#include "Weyl4Scalar.impl.hpp"

#endif /* WEYL4SCALAR_HPP_ */
