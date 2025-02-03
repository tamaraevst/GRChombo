/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEDIAGNOSTICS_HPP_
#define COMPUTEDIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MovingPunctureGauge.hpp"
#include "TensorAlgebra.hpp"
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

template <class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives>
class ComputeDiagnostics : public CCZ4RHS<gauge_t>
{
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    /// CCZ4 variables
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>; 

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = CCZ4Vars::Diff2VarsWithGauge<data_t>;

    struct params_t
    {
        // lapse params:
        double lapse_advec_coeff = 0.; //!< Switches advection terms in
                                       //! the lapse condition on/off
        double lapse_power = 1.; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        double lapse_coeff = 2.; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
        // shift params:
        double shift_Gamma_coeff = 0.75; //!< Gives the F in \f$\partial_t
                                         //!  \beta^i =  F B^i\f$
        double shift_advec_coeff = 0.;   //!< Switches advection terms in the
                                         //! shift condition on/off
        double eta = 1.; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$
    };

     /// Constructor 
    ComputeDiagnostics<gauge_t, deriv_t>(
        params_t a_params, double a_dx,  std::array<double, CH_SPACEDIM> a_center  
    ) : m_params(a_params), m_dx(a_dx), m_center(a_center) {}

    protected:
    params_t m_params;
    double m_dx;
    std::array<double, CH_SPACEDIM> m_center
}
 
template <class gauge_t, class deriv_t>
template <class data_t>
void ComputeDiagnostics<gauge_t, deriv_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell); 
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift); 

    // Spatial 3 metric, g_{ij}
    Tensor<2, data_t> g; 

    FOR2(i,j) g[i][j] = vars.h[i][j] / vars.chi; 

    // Time derivative of the spatial metric, h_{ij} 
    Tensor<2, data_t> ht; 
    data_t divshift = compute_trace(d1.shift);
    auto h_UU = compute_inverse_sym(vars.h);

    FOR(i, j)
    {
        ht[i][j] = advec.h[i][j] - 2.0 * vars.lapse * vars.A[i][j] -
                      (2.0 / GR_SPACEDIM) * vars.h[i][j] * divshift;
        FOR(k)
        {
            ht[i][j] +=
                vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
        }
    }

    // Time derivative of the conformal factor, \chi 
    data_t chit;

    chit = advec.chi +
              (2.0 / GR_SPACEDIM) * vars.chi * (vars.lapse * vars.K - divshift);

    // Time derivative of the spatial metric, g_{ij} 

    Tensor<2, data_t> gt;

    FOR(i,j) 
    {
        gt[i][j] += - 1./(vars.chi * vars.chi) * chit * vars.h[i][j] + 1/(vars.chi) * ht[i][j];
    }

    // d_r derivative of g_{ij}

    Tensor<2, data_t> hr; 
    data_t chir;
    Tensor<2, data_t> gr;

    Coordinates<data_t> coords(current_cell, m_dx, m_center);
    double r = sqrt(coords.x * coords.x + coords.y * coords.y + coords.z * coords.z);

    FOR(i,j)
    {
        hr[i][j] = 1./r * (coord.x * d1.h[0][i][j] + coord.y * d1.h[1][i][j] + coord.z * d1.h[2][i][j]);
    }

    chir = 1./r * (coords.x * d1.chi[0] + coords.y * d1.chi[1] + coords.z * d1.chi[2]);

    FOR(i,j)
    {
        gr[i][j] += - 1./(vars.chi * vars.chi) * chir * vars.h[i][j] + 1/(vars.chi) * hr[i][j];
    }

    // Time derivative of lapse, \alpha

    data_t lapset;

    lapset = m_params.lapse_advec_coeff * advec.lapse -
                    m_params.lapse_coeff *
                        pow(vars.lapse, m_params.lapse_power) *
                        (vars.K - 2 * vars.Theta);
        
    // find d_r derivative of lapse 

    data_t lapser;

    lapser = 1./r * (coords.x * d1.lapse[0] + coords.y * d1.lapse[1] + coords.z * d1.lapse[2]);

    // Time derivative of shift
    Tensor<1, data_t> shiftt;

    FOR(i)
        {
            shiftt[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B[i];
        }

    // d_r derivative of shift
    Tensor<1, data_t> shiftr;

    FOR(i)
    {
        shiftr[i] = 1./r * (coords.x * d1.shift[0][i] + coords.y * d1.shift[1][i] + coords.z * d1.shift[2][i]);
    }

    // DONE WITH ALL VARS

    current_cell.store_vars(g[0][0], c_gxx); // Write the rhs into the output FArrayBox
    current_cell.store_vars(g[0][1], c_gxy); 
    current_cell.store_vars(g[0][2], c_gxz); 
    current_cell.store_vars(g[1][1], c_gyy);
    current_cell.store_vars(g[1][2], c_gyz); 
    current_cell.store_vars(g[2][2], c_gzz); 

    current_cell.store_vars(gt[0][0], c_gtxx);
    current_cell.store_vars(gt[0][1], c_gtxy);
    current_cell.store_vars(gt[0][2], c_gtxz);
    current_cell.store_vars(gt[1][1], c_gtyy);
    current_cell.store_vars(gt[1][2], c_gtyz);
    current_cell.store_vars(gt[2][2], c_gtzz);

    current_cell.store_vars(gr[0][0], c_grxx);
    current_cell.store_vars(gr[0][1], c_grxy);
    current_cell.store_vars(gr[0][2], c_grxz);
    current_cell.store_vars(gr[1][1], c_gryy);
    current_cell.store_vars(gr[1][2], c_gryz);
    current_cell.store_vars(gr[2][2], c_grzz);

    current_cell.store_vars(shiftt[0], c_shifttx);
    current_cell.store_vars(shiftt[1], c_shiftty);
    current_cell.store_vars(shiftt[2], c_shifttz);

    current_cell.store_vars(shiftr[0], c_shiftrx);
    current_cell.store_vars(shiftr[1], c_shiftry);
    current_cell.store_vars(shiftr[2], c_shiftrz);

    current_cell.store_vars(lapset, c_lapset);
    current_cell.store_vars(lapser, c_lapser);
} 

#endif /* COMPUTEDIAGNOSTICS_HPP_ */
