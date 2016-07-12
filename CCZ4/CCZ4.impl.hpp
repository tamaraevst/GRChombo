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
    idx_t<data_t> idx = m_driver.in_idx(ix, iy, iz);

    vars_t<data_t> vars = m_driver.local_vars(idx);

    vars_t<data_t> d1[CH_SPACEDIM];
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        d1[i] = m_deriv.diff1(idx, i);
    }

    vars_t<data_t> d2[CH_SPACEDIM][CH_SPACEDIM];

    // Repeated derivatives
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
       d2[i][i] = m_deriv.diff2(idx, i);
    }

    // Mixed derivatives
    d2[0][1] = m_deriv.mixed_diff2(idx, 1, 0);
    d2[0][2] = m_deriv.mixed_diff2(idx, 2, 0);
    d2[1][2] = m_deriv.mixed_diff2(idx, 2, 1);

    d2[1][0] = d2[0][1];
    d2[2][0] = d2[0][2];
    d2[2][1] = d2[1][2];

    vars_t<data_t> advec = m_deriv.advection(idx, vars.shift);
    vars_t<data_t> dssp = m_deriv.dissipation(idx);

    vars_t<data_t> rhs = rhs_equation(vars, d1, d2, advec);
    
    // TODO: I really do not like this, but cannot think of a better way to do it yet...
    idx_t<data_t> out_idx = m_driver.out_idx(ix, iy, iz);
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_chi])[out_idx]    = rhs.chi      + m_sigma * dssp.chi     ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h11])[out_idx]    = rhs.h[0][0]  + m_sigma * dssp.h[0][0] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h12])[out_idx]    = rhs.h[0][1]  + m_sigma * dssp.h[0][1] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h13])[out_idx]    = rhs.h[0][2]  + m_sigma * dssp.h[0][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h22])[out_idx]    = rhs.h[1][1]  + m_sigma * dssp.h[1][1] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h23])[out_idx]    = rhs.h[1][2]  + m_sigma * dssp.h[1][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_h33])[out_idx]    = rhs.h[2][2]  + m_sigma * dssp.h[2][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_K])[out_idx]      = rhs.K        + m_sigma * dssp.K       ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A11])[out_idx]    = rhs.A[0][0]  + m_sigma * dssp.A[0][0] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A12])[out_idx]    = rhs.A[0][1]  + m_sigma * dssp.A[0][1] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A13])[out_idx]    = rhs.A[0][2]  + m_sigma * dssp.A[0][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A22])[out_idx]    = rhs.A[1][1]  + m_sigma * dssp.A[1][1] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A23])[out_idx]    = rhs.A[1][2]  + m_sigma * dssp.A[1][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_A33])[out_idx]    = rhs.A[2][2]  + m_sigma * dssp.A[2][2] ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Gamma1])[out_idx] = rhs.Gamma[0] + m_sigma * dssp.Gamma[0];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Gamma2])[out_idx] = rhs.Gamma[1] + m_sigma * dssp.Gamma[1];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Gamma3])[out_idx] = rhs.Gamma[2] + m_sigma * dssp.Gamma[2];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Theta])[out_idx]  = rhs.Theta    + m_sigma * dssp.Theta   ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_lapse])[out_idx]  = rhs.lapse    + m_sigma * dssp.lapse   ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_shift1])[out_idx] = rhs.shift[0] + m_sigma * dssp.shift[0];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_shift2])[out_idx] = rhs.shift[1] + m_sigma * dssp.shift[1];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_shift3])[out_idx] = rhs.shift[2] + m_sigma * dssp.shift[2];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_B1])[out_idx]     = rhs.B[0]     + m_sigma * dssp.B[0]    ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_B2])[out_idx]     = rhs.B[1]     + m_sigma * dssp.B[1]    ;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_B3])[out_idx]     = rhs.B[2]     + m_sigma * dssp.B[2]    ;
}

template <class data_t>
auto
CCZ4::rhs_equation(const vars_t<data_t> &vars,
          const vars_t<data_t> (&d1)[CH_SPACEDIM],
          const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
          const vars_t<data_t> &advec
) -> vars_t<data_t>
{
    vars_t<data_t> rhs;

//    Might want to work through the code and eliminate chi divisions where possible to allow chi to go to zero.
//    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto h_UU = TensorAlgebra::compute_inverse(vars.h);
    auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    tensor<1, data_t> Z_over_chi;
    tensor<1, data_t> Z;
    FOR1(i)
    {
       Z_over_chi[i] = 0.5*(vars.Gamma[i] - chris.contracted[i]);
       Z[i] = vars.chi*Z_over_chi[i];
    }

    auto ricci =  CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    data_t divshift = 0;
    data_t Z_dot_d1lapse = 0;
    FOR1(k)
    {
        divshift += d1[k].shift[k];
        Z_dot_d1lapse += Z[k]*d1[k].lapse;
    }

    data_t dlapse_dot_dchi = 0;
    FOR2(m,n)
    {
        dlapse_dot_dchi += h_UU[m][n]*d1[m].lapse*d1[n].chi;
    }

    tensor<2, data_t> covdtilde2lapse;
    tensor<2, data_t> covd2lapse;
    FOR2(k,l)
    {
        covdtilde2lapse[k][l] = d2[k][l].lapse;
        FOR1(m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l]*d1[m].lapse;
        }
        covd2lapse[k][l] = vars.chi*covdtilde2lapse[k][l] + 0.5*(d1[k].lapse*d1[l].chi + d1[k].chi*d1[l].lapse - vars.h[k][l]*dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM/2.0)*dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi*chris.contracted[i]*d1[i].lapse;
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j]*(vars.chi*d2[i][j].lapse + d1[i].lapse*d1[j].chi);
        }
    }


    tensor<2, data_t> A_UU = TensorAlgebra::raise_all(vars.A, h_UU);

    //A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2    = TensorAlgebra::compute_trace(vars.A, A_UU);
    rhs.chi = advec.chi + (2.0/GR_SPACEDIM)*vars.chi*(vars.lapse*vars.K - divshift);
    FOR2(i,j)
    {
        rhs.h[i][j] = advec.h[i][j] - 2.0*vars.lapse*vars.A[i][j] - (2.0/GR_SPACEDIM)*vars.h[i][j]*divshift;
        FOR1(k)
        {
            rhs.h[i][j] += vars.h[k][i]*d1[j].shift[k] + vars.h[k][j]*d1[i].shift[k];
        }
    }

    tensor<2, data_t> Adot_TF;
    FOR2(i,j)
    {
        Adot_TF[i][j] = -covd2lapse[i][j] + vars.chi*vars.lapse*ricci.LL[i][j];
    }
    TensorAlgebra::make_trace_free(Adot_TF, vars.h, h_UU);

    FOR2(i,j)
    {
        rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] + vars.A[i][j]*(vars.lapse*(vars.K - 2*vars.Theta) - (2.0/GR_SPACEDIM)*divshift);
        FOR1(k)
        {
            rhs.A[i][j] += vars.A[k][i]*d1[j].shift[k] + vars.A[k][j]*d1[i].shift[k];
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
            Gammadot[i] += 2*h_UU[i][j]*(vars.lapse*d1[j].Theta - vars.Theta*d1[j].lapse)
                - 2*A_UU[i][j]*d1[j].lapse
                - vars.lapse*((2*(GR_SPACEDIM-1.0)/(double) GR_SPACEDIM)*h_UU[i][j]*d1[j].K + GR_SPACEDIM*A_UU[i][j]*d1[j].chi/vars.chi)
                - (chris.contracted[j] + 2*m_params.kappa3*Z_over_chi[j])*d1[j].shift[i];

            FOR1(k)
            {
                Gammadot[i] += 2*vars.lapse*chris.ULL[i][j][k]*A_UU[j][k]
                    + h_UU[j][k]*d2[j][k].shift[i]
                    + ((GR_SPACEDIM-2.0)/(double) GR_SPACEDIM)*h_UU[i][j]*d2[j][k].shift[k];
            }
        }
    }

    FOR1(i)
    {
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
    }

    const data_t etaDecay = 1;

    rhs.lapse = m_params.lapse_advec_coeff*advec.lapse - 2*vars.lapse*(vars.K - 2*vars.Theta);
    FOR1(i)
    {
        rhs.shift[i] = m_params.shift_advec_coeff*advec.shift[i] + m_params.shift_gamma_coeff*vars.B[i];
        rhs.B[i] = m_params.shift_advec_coeff*advec.B[i] + (1 - m_params.shift_advec_coeff)*advec.Gamma[i] + Gammadot[i] - m_params.beta_driver*etaDecay*vars.B[i];
    }

    return rhs;
}

template <class data_t> template <class arr_t>
CCZ4::vars_t<data_t>::vars_t(const arr_t& in)
{
    chi      = in[c_chi];
    h[0][0]  = in[c_h11];
    h[0][1]  = in[c_h12];
    h[0][2]  = in[c_h13];
    h[1][1]  = in[c_h22];
    h[1][2]  = in[c_h23];
    h[2][2]  = in[c_h33];

    h[1][0]  = h[0][1];
    h[2][0]  = h[0][2];
    h[2][1]  = h[1][2];

    K        = in[c_K];
    A[0][0]  = in[c_A11];
    A[0][1]  = in[c_A12];
    A[0][2]  = in[c_A13];
    A[1][1]  = in[c_A22];
    A[1][2]  = in[c_A23];
    A[2][2]  = in[c_A33];

    A[1][0]  = A[0][1];
    A[2][0]  = A[0][2];
    A[2][1]  = A[1][2];

    Gamma[0] = in[c_Gamma1];
    Gamma[1] = in[c_Gamma2];
    Gamma[2] = in[c_Gamma3];

    Theta    = in[c_Theta];

    lapse    = in[c_lapse];
    shift[0] = in[c_shift1];
    shift[1] = in[c_shift2];
    shift[2] = in[c_shift3];

    B[0]     = in[c_B1];
    B[1]     = in[c_B2];
    B[2]     = in[c_B3];
}

#endif /* CCZ4_IMPL_HPP_ */
