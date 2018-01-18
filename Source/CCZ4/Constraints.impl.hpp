#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

#include "GRInterval.hpp"
#include "MiscUtils.hpp"
#include "VarsTools.hpp"

inline Constraints::Constraints(double dx,
                                double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_cosmological_constant(cosmological_constant)
{
}

template <class data_t>
void Constraints::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Ham, c_Ham);
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1, c_Mom3>());
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
Constraints::constraints_t<data_t> Constraints::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<tensor<1, data_t>> &d1,
    const diff2_vars_t<tensor<2, data_t>> &d2) const
{
    constraints_t<data_t> out;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

    auto A_UU = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

    out.Ham = ricci.scalar +
              (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
    out.Ham -= 2 * m_cosmological_constant;

    tensor<2, data_t> covd_A[CH_SPACEDIM];
    FOR3(i, j, k)
    {
        covd_A[i][j][k] = d1.A[j][k][i];
        FOR1(l)
        {
            covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                               chris.ULL[l][i][k] * vars.A[l][j];
        }
    }

    FOR1(i) { out.Mom[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM; }
    FOR3(i, j, k)
    {
        out.Mom[i] += h_UU[j][k] *
                      (covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] *
                                             d1.chi[k] / (2 * chi_regularised));
    }

    return out;
}

template <class data_t>
template <typename mapping_function_t>
void Constraints::Vars<data_t>::enum_mapping(
    mapping_function_t mapping_function)
{
    using namespace VarsTools; // define_enum_mapping is part of VarsTools
    define_enum_mapping(mapping_function, c_chi, chi);
    define_enum_mapping(mapping_function, c_K, K);
    define_enum_mapping(mapping_function, GRInterval<c_Gamma1, c_Gamma3>(),
                        Gamma);
    define_symmetric_enum_mapping(mapping_function, GRInterval<c_h11, c_h33>(),
                                  h);
    define_symmetric_enum_mapping(mapping_function, GRInterval<c_A11, c_A33>(),
                                  A);
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
