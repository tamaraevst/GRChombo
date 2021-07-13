/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SCALARFIELDMODIFIEDGR_HPP_)
#error "This file should only be included through ScalarFieldModifiedGR.hpp"
#endif

#ifndef SCALARFIELDMODIFIEDGR_IMPL_HPP_
#define SCALARFIELDMODIFIEDGR_IMPL_HPP_

// the RHS excluding the potential terms
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void matter_rhs_modified_gr(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    auto starR_R = CCZ4GeometryModifieddGR::compute_chern_simons(vars, d1, d2, h_UU, chris);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi = vars.lapse * vars.Pi + advec.phi - starR_R;
    rhs.Pi = vars.lapse * vars.K * vars.Pi + advec.Pi;

    FOR2(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi += h_UU[i][j] * (-0.5 * d1.chi[j] * vars.lapse * d1.phi[i] +
                                vars.chi * vars.lapse * d2.phi[i][j] +
                                vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR1(k)
        {
            rhs.Pi += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi[k];
        }
    }
}

#endif /* SCALARFIELD_IMPL_HPP_ */
