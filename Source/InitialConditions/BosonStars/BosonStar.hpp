/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTAR_HPP_
#define BOSONSTAR_HPP_

// TODO: uncomment these includes when incorporating into GRChombo.
//#include "Cell.hpp"
//#include "Coordinates.hpp"
//#include "MatterCCZ4.hpp"
//#include "ComplexScalarField.hpp"
//#include "Tensor.hpp"
//#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
//#include "VarsTools.hpp"
//#include "simd.hpp"
#include "ComplexPotential.hpp"

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BosonStar
{
  public:
    //! A structure for the input params for the boson star
    struct params_t
    {
        double central_amplitude_CSF; //!< Central amplitude of the star
        // TODO: uncomment when incorporating into GRCHombo
        //std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
        double abs_error; //!< absolute error tolerance for the 1D ODE integrator
        double rel_error; //!< relative error tolerance for the 1D ODE integrator
        double initial_step_size; //!< initial step size for the 1D ODE integrator
        double max_radius; /*!< the maximum (rescaled) radius the 1D ODE
        integrator will attempt to integrate up to. This should be set higher
        than it will actually be able to reach.*/
        double binary_search_tol; /*!< This is the tolerance to which the shooting
        parameter (the central value of alpha) is found. */
        int max_binary_search_iter = 1e3; //!< the maximum number of binary search iterations.
    };

    //! The constructor
    BosonStar(params_t a_params_CSF,
        Potential::params_t a_params_potential, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    //<-(remove)template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    const params_t m_params_CSF; //!< The complex scalar field params
    const Potential::params_t m_params_potential; //!< The potential params

};

#include "BosonStar.impl.hpp"

#endif /* BOSONSTAR_HPP_ */
