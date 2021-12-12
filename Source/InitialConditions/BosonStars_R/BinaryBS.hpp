/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBS_HPP_
#define BINARYBS_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"
#include "BosonStarSolution.hpp"
#include <vector>
#include "parstream.H" //gives pout
#include "WeightFunction.hpp"

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BinaryBS
{

public:
    //! The constructor
    // BosonStar(BosonStar_params_t a_params_BosonStar, BosonStar_params_t a_params_BosonStar2,
    //     Potential::params_t a_params_potential, double a_G_Newton, double a_dx,
    //     int a_verbosity);

    BinaryBS(BosonStar_params_t a_bosonstar_params, BosonStar_params_t a_bosonstar2_params,
        Potential::params_t a_params_potential, double a_G_Newton, double a_dx,
        int a_verbosity);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

    BosonStarSolution m_bosonstar;
    BosonStarSolution m_bosonstar2; /*<
    The object that stores the solution found by the 1d ODE integrator */


protected:
    double m_dx;
    double m_G_Newton;
    // BosonStar_params_t m_params_BosonStar; //!< The complex scalar field params
    // BosonStar_params_t m_params_BosonStar2;
    Potential::params_t m_params_potential; //!< The potential params
    int m_verbosity;

    BosonStar_params_t m_bosonstar_params;
    BosonStar_params_t m_bosonstar2_params;

};

#include "BinaryBS.impl.hpp"

#endif /* BINARYBS_HPP_ */
