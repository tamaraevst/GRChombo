/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPUTEMODIFIEDSCALARS_HPP
#define COMPUTEMODIFIEDSCALARS_HPP

#include "Cell.hpp"
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "BSSNVars.hpp"
#include "Tensor.hpp"
#include "ModifiedScalarsGeometry.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include "DebuggingTools.hpp"

//! This class computes Chern Simons and Gauss Bonnet scalars
class ComputeModifiedScalars
{
  public:
    /// BSSN variables
    template <class data_t> using Vars = BSSNVars::VarsWithGauge<data_t>;

    /// BSSN variables
     template <class data_t>
    using Diff2Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;
    // template <class data_t>
    // using Diff2Vars = BSSNVars::Diff2VarsWithGauge<data_t>;

    //! Constructor
    ComputeModifiedScalars(const std::array<double, CH_SPACEDIM> a_center,
                     const double a_dx,
                     double a_gamma_amplitude, 
                     double a_beta_amplitude)
                     : m_center(a_center),
                     m_dx(a_dx),
                     m_gamma_amplitude(a_gamma_amplitude), 
                     m_beta_amplitude(a_beta_amplitude),
                     m_deriv(a_dx){}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);

        // Calculate modified scalars
        ModifiedScalarsGeometry mod_geom;
        auto mod_scalars = mod_geom.mod_scalars(vars, d1, d2, h_UU, chris);

        current_cell.store_vars(mod_scalars.GB_scalar, c_gaussbonnet);
        current_cell.store_vars(mod_scalars.CS_scalar, c_chernsimons);
    }

  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    double m_beta_amplitude;
    double m_gamma_amplitude;
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
};


#endif /* COMPUTEMODIFIEDSCALARS_HPP */