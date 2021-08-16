/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPUTEMODIFIEDSCALARS_HPP)
#error "This file should only be included through ComputeModifiedScalars.hpp"
#endif

#ifndef COMPUTEMODIFIEDSCALARS_IMPL_HPP
#define COMPUTEMODIFIEDSCALARS_IMPL_HPP

ComputeModifiedScalars::ComputeModifiedScalars(
    const std::array<double, CH_SPACEDIM> &a_center,
    const double a_dx, const int a_var_enum)
    : m_center(a_center), m_dx(a_dx), m_deriv(a_dx),
      m_var_enum(a_var_enum)
{
}

template <class data_t>
void ComputeModifiedScalars::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Calculate EM tensor
    auto mod_scalar = CCZ4GeometryModifiedGR::compute_modified_scalars(vars, d1, d2, h_UU, chris);

    current_cell.store_vars(mod_scalar.starR_R, m_var_enum);
}

#endif /* COMPUTEMODIFIEDSCALARS_IMPL_HPP_ */