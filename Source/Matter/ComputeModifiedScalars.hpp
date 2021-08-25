/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEMODIFIEDSCALARS_HPP
#define COMPUTEMODIFIEDSCALARS_HPP

#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "BSSNVars.hpp"
#include "Tensor.hpp"
#include "CCZ4GeometryModifiedGR.hpp"
#include "simd.hpp"
#include "Coordinates.hpp"
#include "DebuggingTools.hpp"

//! This class computes Chern Simons and Gauss Bonnet scalars
class ComputeModifiedScalars
{
  public:
    /// CCZ4 variables
    template <class data_t> using MetricVars = BSSNVars::VarsWithGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsWithGauge<data_t>;

    /// Vars object for Constraints
    template <class data_t> struct Vars
    {
        data_t RGB;
        data_t starR_R;

        // template <typename mapping_function_t>
        // void enum_mapping(mapping_function_t mapping_function)
        // {
        //     using namespace VarsTools;
        //     define_enum_mapping(mapping_function, c_chernsimons, starR_R);
        //     define_enum_mapping(mapping_function, c_gaussbonnet, RGB);
        // }
    };

    //! Constructor
    ComputeModifiedScalars(const std::array<double, CH_SPACEDIM> a_center,
                     const double a_dx,
                     double a_gamma_amplitude, 
                     double a_beta_amplitude, 
                     const int a_c_chernsimons, 
                     const int a_c_gaussbonnet)
                     : m_center(a_center),
                     m_dx(a_dx),
                     m_gamma_amplitude(a_gamma_amplitude), 
                     m_beta_amplitude(a_beta_amplitude), 
                     m_c_chernsimons(a_c_chernsimons), 
                     m_c_gaussbonnet(a_c_gaussbonnet),
                     m_deriv(a_dx){}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<MetricVars>();
        const auto d1 = m_deriv.template diff1<MetricVars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);
        
        // Calculate modified scalars
        Vars<data_t> out = modified_scalars(vars, d1, d2, h_UU);

        // Get the coordinates  
        const Coordinates<double> coords(current_cell, m_dx, m_center);

        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        
        if (x>6.0 || y>6.0 ||z>6.0)
        {
          out.starR_R = 0.0;
        }

        DEBUG_OUT(x);
        DEBUG_OUT(y);
        DEBUG_OUT(z);
        DEBUG_OUT(out.starR_R);

        store_vars(out, current_cell);
    }

  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_c_chernsimons;
    const int m_c_gaussbonnet;
    double m_beta_amplitude;
    double m_gamma_amplitude;

    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

    template <class data_t>
    void store_vars(Vars<data_t> &out, Cell<data_t> &current_cell) const
    { 
      current_cell.store_vars(out.RGB, m_c_gaussbonnet);
      current_cell.store_vars(out.starR_R, m_c_chernsimons);
    }
  
    // Function for computing Gauss Bonnet and Chern Simons scalars!
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    Vars<data_t> modified_scalars(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU) const
    {
        using namespace TensorAlgebra;
        Vars<data_t> out;

        CCZ4GeometryModifiedGR ccz4mod;

        // const data_t x = coords.x;
        // const double y = coords.y;
        // const double z = coords.z;

        const auto E_ij = ccz4mod.compute_chern_simons_electric_term(vars, d1, d2, h_UU);
        const auto B_ij = ccz4mod.compute_magnetic_term(vars, d1, d2, h_UU);

        FOR4(i, j, k, l)
        {
            out.starR_R = - 8.0 * vars.chi * vars.chi * h_UU[k][i] * h_UU[l][j] * B_ij[k][l] * E_ij[i][j];
        }

        out.RGB = ccz4mod.GB_scalar(vars, d1, d2, h_UU);
        // if (y>10.0 || z>10.0)
        // {   
        //     out.starR_R = 0.0;
        // }
        // else{
        //     //Finally compute *RR
        //     FOR4(i, j, k, l)
        // {
        //     out.starR_R = - 8.0 * vars.chi * vars.chi * h_UU[k][i] * h_UU[l][j] * B_ij[k][l] * E_ij[i][j];
        // }

        // out.RGB = ccz4mod.GB_scalar(vars, d1, d2, h_UU);
        // }
        DEBUG_OUT(out.starR_R);
        return out;
    }
};


#endif /* COMPUTEMODIFIEDSCALARS_HPP */