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
#include "CCZ4GeometryModifiedGR.hpp"
#include "simd.hpp"

class ComputeModifiedScalars
{
  public:

    //! Constructor
    ComputeModifiedScalars(const std::array<double, CH_SPACEDIM> &a_center,
                     const double a_dx, const int a_var_enum);

    template <class data_t> void compute(Cell<data_t> current_cell) const;

    const std::array<double, CH_SPACEDIM> m_center;
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_var_enum;
};

#include "ComputeModifiedScalars.impl.hpp"

#endif /* COMPUTEMODIFIEDSCALARS_HPP */