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
#include "MatterCCZ4.hpp"
#include "CCZ4GeometryModifiedGR.hpp"
#include "simd.hpp"

//! This class computes Chern Simons and Gauss Bonnet scalars
template <class matter_t> class ComputeModifiedScalars
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    ComputeModifiedScalars(matter_t a_matter,
                     const double a_dx,
                     double a_gamma_amplitude, 
                     double a_beta_amplitude, 
                     const int a_c_chernsimons, 
                     const int a_c_gaussbonnet)
                     : m_matter(a_matter), 
                     m_dx(a_dx),
                     m_gamma_amplitude(a_gamma_amplitude), 
                     m_beta_amplitude(a_beta_amplitude), 
                     m_c_chernsimons(a_c_chernsimons), 
                     m_c_gaussbonnet(a_c_gaussbonnet),
                     m_deriv(a_dx){}

  protected:
    matter_t m_matter;
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_c_chernsimons;
    const int m_c_gaussbonnet;
    double m_beta_amplitude;
    double m_gamma_amplitude;

  public:
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto matter_vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(matter_vars.h);
        const auto chris = compute_christoffel(d1.h, h_UU);

        CCZ4GeometryModifiedGR ccz4mod(m_gamma_amplitude, m_beta_amplitude);
        
        // Calculate modified scalars
        const auto modified_scalars = ccz4mod.add_modified_scalars(matter_vars, d1, d2, h_UU, chris);

        current_cell.store_vars(modified_scalars.RGB, m_c_gaussbonnet);
        current_cell.store_vars(modified_scalars.starR_R, m_c_chernsimons);
    }
};


#endif /* COMPUTEMODIFIEDSCALARS_HPP */