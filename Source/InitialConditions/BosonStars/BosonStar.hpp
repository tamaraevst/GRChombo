/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTAR_HPP_
#define BOSONSTAR_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "ComplexPotential.hpp"
#include <boost/numerics/odeint.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BosonStar
{
  public:
    //! A structure for the input params for the boson star
    struct params_t
    {
        double central_amplitude_CSF; //!< Central amplitude of the star
        std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
    };

    //! The constructor
    BosonStar(params_t a_params_CSF,
        ComplexPotential::params_t a_params_potential, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    const params_t m_params_CSF; //!< The complex scalar field params
    const ComplexPotential::params_t m_params_potential; //!< The potential params

};

#include "BosonStar.impl.hpp"

#endif /* BOSONSTAR_HPP_ */
