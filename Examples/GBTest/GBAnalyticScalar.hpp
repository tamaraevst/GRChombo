/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GBANALYTICSCALAR_HPP_
#define GBANALYTICSCALAR_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

class GBAnalyticScalar
{
    // Use the variable definition in CCZ4
    template <class data_t>
    using Vars = ADMConformalVars::VarsWithGauge<data_t>;

    protected:
        const double m_dx;
        const std::array<double, CH_SPACEDIM> m_center; //!< The center of the Kerr BH

    public:
        GBAnalyticScalar(const double a_dx, const std::array<double, CH_SPACEDIM> a_center) : m_dx(a_dx), m_center(a_center){};
    
        template <class data_t> void compute(Cell<data_t> current_cell) const   
        {
            // The cartesian variables and coords
            Vars<data_t> vars;
            Coordinates<data_t> coords(current_cell, m_dx, m_center);

            //This is just for the computation of the analytic phi expression on Schwarzschild background. We fix M and beta.
            double M = 1.0;
            double beta = 0.5;
        
            data_t r =  sqrt(vars.chi) * coords.get_radius();
            data_t xx = pow((1.0 + M / (2.0 * r)), 2.0) * r / M;
            data_t phi_analytic = ((2.0 * beta) / (M * M)) * (1.0 / xx + 1.0 / (xx * xx) + (4.0 / 3.0) * 1.0 / (xx * xx * xx));

            current_cell.store_vars(phi_analytic, c_phianalytic);            

        }
};

#endif /* GBANALYTICSCALAR_HPP_ */
