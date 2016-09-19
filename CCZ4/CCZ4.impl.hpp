#if !defined(CCZ4_HPP_)
#error "This file should only be included through CCZ4.hpp"
#endif

#ifndef CCZ4_IMPL_HPP_
#define CCZ4_IMPL_HPP_

#define COVARIANTZ4

CCZ4::CCZ4(const FABDriverBase& driver, params_t params, double dx, double sigma, double cosmological_constant) :
    m_params (params),
    m_sigma (sigma),
    m_cosmological_constant (cosmological_constant),
    m_driver (driver),
    m_deriv (dx, m_driver)
{}

template <class data_t>
void
CCZ4::compute(int ix, int iy, int iz)
{
    idx_t<data_t> idx = m_driver.in_idx(ix, iy, iz); //The current location in the flattened input FArrayBox

    vars_t<data_t> vars;
    m_driver.local_vars(vars,idx);

    vars_t< tensor<1, data_t> > d1;
    FOR1(idir) m_deriv.diff1(d1, idx, idir);

    vars_t< tensor<2,data_t> > d2;
    // Repeated derivatives
    FOR1(idir) m_deriv.diff2(d2, idx, idir);
    // Mixed derivatives
    // Note: no need to symmetrise explicitely, this is done in mixed_diff2
    m_deriv.mixed_diff2(d2, idx, 1, 0);
    m_deriv.mixed_diff2(d2, idx, 2, 0);
    m_deriv.mixed_diff2(d2, idx, 2, 1);

    vars_t<data_t> advec;
    advec.assign(0.);
    FOR1(idir) m_deriv.add_advection(advec, idx, vars.shift[idir], idir);

    vars_t<data_t> rhs = rhs_equation(vars, d1, d2, advec);

    FOR1(idir) m_deriv.add_dissipation(rhs, m_sigma, idx,idir);

    //Write the rhs into the output FArrayBox
    idx_t<data_t> out_idx = m_driver.out_idx(ix, iy, iz); //The current location in the flattened output FArraBox
    m_driver.store_vars(rhs, out_idx);
}

template <class data_t>
auto
CCZ4::rhs_equation(const vars_t<data_t> &vars,
          const vars_t< tensor<1,data_t> >& d1,
          const vars_t< tensor<2,data_t> >& d2,
          const vars_t<data_t> &advec
) -> vars_t<data_t>
{
    vars_t<data_t> rhs;

//    Might want to work through the code and eliminate chi divisions where possible to allow chi to go to zero.
//    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    using namespace TensorAlgebra;

    auto h_UU = compute_inverse(vars.h);
    auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    tensor<1, data_t> Z_over_chi;
    tensor<1, data_t> Z;
    FOR1(i)
    {
       Z_over_chi[i] = 0.5*(vars.Gamma[i] - chris.contracted[i]);
       Z[i] = vars.chi*Z_over_chi[i];
    }

    auto ricci =  CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    data_t divshift = compute_trace(d1.shift);
    data_t Z_dot_d1lapse = compute_dot_product(Z,d1.lapse);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse,d1.chi,h_UU);

    tensor<2, data_t> covdtilde2lapse;
    tensor<2, data_t> covd2lapse;
    FOR2(k,l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR1(m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l]*d1.lapse[m];
        }
        covd2lapse[k][l] = vars.chi*covdtilde2lapse[k][l] + 0.5*(d1.lapse[k]*d1.chi[l] + d1.chi[k]*d1.lapse[l] - vars.h[k][l]*dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM/2.0)*dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi*chris.contracted[i]*d1.lapse[i];
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j]*(vars.chi*d2.lapse[i][j] + d1.lapse[i]*d1.chi[j]);
        }
    }


    tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    //A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2    = compute_trace(vars.A, A_UU);
    rhs.chi = advec.chi + (2.0/GR_SPACEDIM)*vars.chi*(vars.lapse*vars.K - divshift);
    FOR2(i,j)
    {
        rhs.h[i][j] = advec.h[i][j] - 2.0*vars.lapse*vars.A[i][j] - (2.0/GR_SPACEDIM)*vars.h[i][j]*divshift;
        FOR1(k)
        {
            rhs.h[i][j] += vars.h[k][i]*d1.shift[k][j] + vars.h[k][j]*d1.shift[k][i];
        }
    }

    tensor<2, data_t> Adot_TF;
    FOR2(i,j)
    {
        Adot_TF[i][j] = -covd2lapse[i][j] + vars.chi*vars.lapse*ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, vars.h, h_UU);

    FOR2(i,j)
    {
        rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] + vars.A[i][j]*(vars.lapse*(vars.K - 2*vars.Theta) - (2.0/GR_SPACEDIM)*divshift);
        FOR1(k)
        {
            rhs.A[i][j] += vars.A[k][i]*d1.shift[k][j] + vars.A[k][j]*d1.shift[k][i];
            FOR1(l)
            {
                rhs.A[i][j] -= 2*vars.lapse*h_UU[k][l]*vars.A[i][k]*vars.A[l][j];
            }
        }
    }

#ifdef COVARIANTZ4
    data_t kappa1_lapse = m_params.kappa1;
#else
    data_t kappa1_lapse = m_params.kappa1 * vars.lapse;
#endif

    rhs.Theta = advec.Theta + 0.5*vars.lapse*(ricci.scalar - tr_A2 + ((GR_SPACEDIM-1.0)/(double) GR_SPACEDIM)*vars.K*vars.K - 2*vars.Theta*vars.K) - 0.5*vars.Theta*kappa1_lapse*((GR_SPACEDIM+1) + m_params.kappa2*(GR_SPACEDIM-1)) - Z_dot_d1lapse;

    rhs.Theta += - vars.lapse * m_cosmological_constant;

    rhs.K = advec.K + vars.lapse*(ricci.scalar + vars.K*(vars.K - 2*vars.Theta) ) - kappa1_lapse*GR_SPACEDIM*(1+m_params.kappa2)*vars.Theta - tr_covd2lapse;

    tensor<1, data_t> Gammadot;
    FOR1(i)
    {
        Gammadot[i] = (2.0/GR_SPACEDIM)*(divshift*(chris.contracted[i] + 2*m_params.kappa3*Z_over_chi[i]) - 2*vars.lapse*vars.K*Z_over_chi[i]) - 2*kappa1_lapse*Z_over_chi[i];
        FOR1(j)
        {
            Gammadot[i] += 2*h_UU[i][j]*(vars.lapse*d1.Theta[j] - vars.Theta*d1.lapse[j])
                - 2*A_UU[i][j]*d1.lapse[j]
                - vars.lapse*((2*(GR_SPACEDIM-1.0)/(double) GR_SPACEDIM)*h_UU[i][j]*d1.K[j] + GR_SPACEDIM*A_UU[i][j]*d1.chi[j]/vars.chi)
                - (chris.contracted[j] + 2*m_params.kappa3*Z_over_chi[j])*d1.shift[i][j];

            FOR1(k)
            {
                Gammadot[i] += 2*vars.lapse*chris.ULL[i][j][k]*A_UU[j][k]
                    + h_UU[j][k]*d2.shift[i][j][k]
                    + ((GR_SPACEDIM-2.0)/(double) GR_SPACEDIM)*h_UU[i][j]*d2.shift[k][j][k];
            }
        }
    }

    FOR1(i)
    {
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
    }

    const data_t etaDecay = 1.;

    rhs.lapse = m_params.lapse_advec_coeff*advec.lapse - 2*vars.lapse*(vars.K - 2*vars.Theta);
    FOR1(i)
    {
        rhs.shift[i] = m_params.shift_advec_coeff*advec.shift[i] + m_params.shift_gamma_coeff*vars.B[i];
        rhs.B[i] = m_params.shift_advec_coeff*advec.B[i] + (1 - m_params.shift_advec_coeff)*advec.Gamma[i] + Gammadot[i] - m_params.beta_driver*etaDecay*vars.B[i];
    }

    return rhs;
}

template <class data_t>
CCZ4::vars_t<data_t>::vars_t()
{
    m_assignment_ptrs[c_chi].push_back(&chi);

    m_assignment_ptrs[c_h11].push_back(&h[0][0]);
    m_assignment_ptrs[c_h12].push_back(&h[0][1]);
    m_assignment_ptrs[c_h12].push_back(&h[1][0]);
    m_assignment_ptrs[c_h13].push_back(&h[0][2]);
    m_assignment_ptrs[c_h13].push_back(&h[2][0]);
    m_assignment_ptrs[c_h22].push_back(&h[1][1]);
    m_assignment_ptrs[c_h23].push_back(&h[1][2]);
    m_assignment_ptrs[c_h23].push_back(&h[2][1]);
    m_assignment_ptrs[c_h33].push_back(&h[2][2]);

    m_assignment_ptrs[c_K].push_back(&K);

    m_assignment_ptrs[c_A11].push_back(&A[0][0]);
    m_assignment_ptrs[c_A12].push_back(&A[0][1]);
    m_assignment_ptrs[c_A12].push_back(&A[1][0]);
    m_assignment_ptrs[c_A13].push_back(&A[0][2]);
    m_assignment_ptrs[c_A13].push_back(&A[2][0]);
    m_assignment_ptrs[c_A22].push_back(&A[1][1]);
    m_assignment_ptrs[c_A23].push_back(&A[1][2]);
    m_assignment_ptrs[c_A23].push_back(&A[2][1]);
    m_assignment_ptrs[c_A33].push_back(&A[2][2]);

    m_assignment_ptrs[c_Gamma1].push_back(&Gamma[0]);
    m_assignment_ptrs[c_Gamma2].push_back(&Gamma[1]);
    m_assignment_ptrs[c_Gamma3].push_back(&Gamma[2]);

    m_assignment_ptrs[c_Theta].push_back(&Theta);

    m_assignment_ptrs[c_lapse].push_back(&lapse);
    m_assignment_ptrs[c_shift1].push_back(&shift[0]);
    m_assignment_ptrs[c_shift2].push_back(&shift[1]);
    m_assignment_ptrs[c_shift3].push_back(&shift[2]);

    m_assignment_ptrs[c_B1].push_back(&B[0]);
    m_assignment_ptrs[c_B2].push_back(&B[1]);
    m_assignment_ptrs[c_B3].push_back(&B[2]);
}

#endif /* CCZ4_IMPL_HPP_ */
