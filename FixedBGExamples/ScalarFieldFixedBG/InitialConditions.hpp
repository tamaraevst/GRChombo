/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_amplitude;
    const std::array<double, CH_SPACEDIM> m_center;
    const IsotropicKerrFixedBG::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const double a_amplitude,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const IsotropicKerrFixedBG::params_t a_bg_params,
                      const double a_dx)
        : m_dx(a_dx), m_amplitude(a_amplitude), m_center(a_center), 
          m_bg_params(a_bg_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        IsotropicKerrFixedBG boosted_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        boosted_bh.compute_metric_background(metric_vars, current_cell);
        const data_t det_gamma =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        data_t phi = m_amplitude;
        
        // Store the initial values of the variables
        current_cell.store_vars(phi, c_phi);
        current_cell.store_vars(0.0, c_Pi);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
