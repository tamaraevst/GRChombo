#include "CCZ4.hpp"

#define COVARIANTZ4

CCZ4::CCZ4(params_t params, double dx, double sigma) :
    m_params (params),
    m_dx (dx),
    m_sigma (sigma)
{}

void
CCZ4::execute(const FArrayBox& in, FArrayBox& out)
{
    // dataPtr in Chombo does CH_assert bound check
    // which we don't want to do in a loop
    for (int i = 0; i < c_NUM; ++i)
    {
        m_in_ptr[i] = in.dataPtr(i);
        m_out_ptr[i] = out.dataPtr(i);
    }

    m_in_lo = in.loVect();
    m_in_hi = in.hiVect();
    m_stride[0] = 1;
    m_stride[1] = m_in_hi[0]-m_in_lo[0]+1;
    m_stride[2] = (m_in_hi[1]-m_in_lo[1]+1)*m_stride[1];

    m_out_lo = out.loVect();
    m_out_hi = out.hiVect(); 
    m_out_stride[0] = 1;
    m_out_stride[1] = m_out_hi[0]-m_out_lo[0]+1;
    m_out_stride[2] = (m_out_hi[1]-m_out_lo[1]+1)*m_out_stride[1];

#pragma omp parallel for default(shared) collapse(2)
    for (int z = m_out_lo[2]; z <= m_out_hi[2]; ++z)
    for (int y = m_out_lo[1]; y <= m_out_hi[1]; ++y)
    {
        int x_simd_max = m_out_lo[0] + simd<double>::simd_len * (((m_out_hi[0] - m_out_lo[0] + 1) / simd<double>::simd_len) - 1);

        // SIMD LOOP
        #pragma novector
        for (int x = m_out_lo[0]; x <= x_simd_max; x += simd<double>::simd_len)
        {
            compute<simd<double> >(x, y, z);
        }

        // REMAINDER LOOP
        #pragma novector
        for (int x = x_simd_max + simd<double>::simd_len; x <= m_out_hi[0]; ++x)
        {
            compute<double>(x, y, z);
        }
    }

}

template <class data_t>
void
CCZ4::compute(int x, int y, int z)
{
    const int idx = m_stride[2]*(z-m_in_lo[2]) + m_stride[1]*(y-m_in_lo[1]) + (x-m_in_lo[0]);

    vars_t<data_t> vars;
    local_vars(idx, vars);
    
    vars_t<data_t> d1[CH_SPACEDIM];
    for (int i = 0; i < 3; ++i)
    {
        diff1(idx, m_stride[i], d1[i]);
    }

    vars_t<data_t> d2[CH_SPACEDIM][CH_SPACEDIM];

    // Repeated derivatives
    for (int i = 0; i < 3; ++i)
    {
        diff2(idx, m_stride[i], d2[i][i]);
    }

    // Mixed derivatives
    mixed_diff2(idx, m_stride[1], m_stride[0], d2[0][1]);
    mixed_diff2(idx, m_stride[2], m_stride[0], d2[0][2]);
    mixed_diff2(idx, m_stride[2], m_stride[1], d2[1][2]);

    d2[1][0] = d2[0][1];
    d2[2][0] = d2[0][2];
    d2[2][1] = d2[1][2];

    vars_t<data_t> advec;
    advection(idx, vars.shift, advec);

    vars_t<data_t> dssp;
    dissipation(idx, dssp);

    vars_t<data_t> rhs = {};
    rhs_equation(vars, d1, d2, advec, rhs);

    // TODO: I really do not like this, but cannot think of a better way to do it yet...
    const int out_idx = m_out_stride[2]*(z-m_out_lo[2]) + m_out_stride[1]*(y-m_out_lo[1]) + (x-m_out_lo[0]);
    SIMDIFY<data_t>(m_out_ptr[c_chi])[out_idx]    = rhs.chi      + m_sigma * dssp.chi;
    SIMDIFY<data_t>(m_out_ptr[c_h11])[out_idx]    = rhs.h[0][0]  + m_sigma * dssp.h[0][0];
    SIMDIFY<data_t>(m_out_ptr[c_h12])[out_idx]    = rhs.h[0][1]  + m_sigma * dssp.h[0][1];
    SIMDIFY<data_t>(m_out_ptr[c_h13])[out_idx]    = rhs.h[0][2]  + m_sigma * dssp.h[0][2];
    SIMDIFY<data_t>(m_out_ptr[c_h22])[out_idx]    = rhs.h[1][1]  + m_sigma * dssp.h[1][1];
    SIMDIFY<data_t>(m_out_ptr[c_h23])[out_idx]    = rhs.h[1][2]  + m_sigma * dssp.h[1][2];
    SIMDIFY<data_t>(m_out_ptr[c_h33])[out_idx]    = rhs.h[2][2]  + m_sigma * dssp.h[2][2];
    SIMDIFY<data_t>(m_out_ptr[c_K])[out_idx]      = rhs.K        + m_sigma * dssp.K;
    SIMDIFY<data_t>(m_out_ptr[c_A11])[out_idx]    = rhs.A[0][0]  + m_sigma * dssp.A[0][0];
    SIMDIFY<data_t>(m_out_ptr[c_A12])[out_idx]    = rhs.A[0][1]  + m_sigma * dssp.A[0][1];
    SIMDIFY<data_t>(m_out_ptr[c_A13])[out_idx]    = rhs.A[0][2]  + m_sigma * dssp.A[0][2];
    SIMDIFY<data_t>(m_out_ptr[c_A22])[out_idx]    = rhs.A[1][1]  + m_sigma * dssp.A[1][1];
    SIMDIFY<data_t>(m_out_ptr[c_A23])[out_idx]    = rhs.A[1][2]  + m_sigma * dssp.A[1][2];
    SIMDIFY<data_t>(m_out_ptr[c_A33])[out_idx]    = rhs.A[2][2]  + m_sigma * dssp.A[2][2];
    SIMDIFY<data_t>(m_out_ptr[c_Gamma1])[out_idx] = rhs.Gamma[0] + m_sigma * dssp.Gamma[0];
    SIMDIFY<data_t>(m_out_ptr[c_Gamma2])[out_idx] = rhs.Gamma[1] + m_sigma * dssp.Gamma[1];
    SIMDIFY<data_t>(m_out_ptr[c_Gamma3])[out_idx] = rhs.Gamma[2] + m_sigma * dssp.Gamma[2];
    SIMDIFY<data_t>(m_out_ptr[c_Theta])[out_idx]  = rhs.Theta    + m_sigma * dssp.Theta;
    SIMDIFY<data_t>(m_out_ptr[c_lapse])[out_idx]  = rhs.lapse    + m_sigma * dssp.lapse;
    SIMDIFY<data_t>(m_out_ptr[c_shift1])[out_idx] = rhs.shift[0] + m_sigma * dssp.shift[0];
    SIMDIFY<data_t>(m_out_ptr[c_shift2])[out_idx] = rhs.shift[1] + m_sigma * dssp.shift[1];
    SIMDIFY<data_t>(m_out_ptr[c_shift3])[out_idx] = rhs.shift[2] + m_sigma * dssp.shift[2];
    SIMDIFY<data_t>(m_out_ptr[c_B1])[out_idx]     = rhs.B[0]     + m_sigma * dssp.B[0];
    SIMDIFY<data_t>(m_out_ptr[c_B2])[out_idx]     = rhs.B[1]     + m_sigma * dssp.B[1];
    SIMDIFY<data_t>(m_out_ptr[c_B3])[out_idx]     = rhs.B[2]     + m_sigma * dssp.B[2];
}

template <class data_t>
void
CCZ4::rhs_equation(const vars_t<data_t> &vars,
          const vars_t<data_t> (&d1)[CH_SPACEDIM],
          const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
          const vars_t<data_t> &advec,
          vars_t<data_t> &rhs
)
{
    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    auto inv = compute_inverse_metric(vars);
    chris_t chris(vars, d1, inv);

    tensor<1, data_t> Z_over_chi;
    tensor<1, data_t> Z;
    {
        FOR1(i)
        {
            Z_over_chi[i] = 0.5*(vars.Gamma[i] - chrisvec[i]);
            Z[i] = vars.chi*Z_over_chi[i];
        }
    }

    ricciZ_t ricci(vars, d2, d1, chris, h_UU, Z_over_chi);


    data_t divshift = 0;
    data_t Z_dot_d1lapse = 0;
    {
        FOR1(k)
        {
            divshift += d1[k].shift[k]; 
            Z_dot_d1lapse += Z[k]*d1[k].lapse;
        }
    }

    data_t dchi_dot_dchi = 0;
    data_t dlapse_dot_dchi = 0;
    {
        FOR2(m,n)
        {
            dchi_dot_dchi += h_UU[m][n]*d1[m].chi*d1[n].chi;
            dlapse_dot_dchi += h_UU[m][n]*d1[m].lapse*d1[n].chi;
        }
    }

    tensor<2, data_t> covdtilde2lapse;
    tensor<2, data_t> covd2lapse;
    {
        FOR2(k,l)
        {
            covdtilde2chi[k][l] = d2[k][l].chi;
            covdtilde2lapse[k][l] = d2[k][l].lapse;
            FOR1(m)
            {
                covdtilde2chi[k][l] -= chris[m][k][l]*d1[m].chi;
                covdtilde2lapse[k][l] -= chris[m][k][l]*d1[m].lapse;
            }
            covd2lapse[k][l] = vars.chi*covdtilde2lapse[k][l] + 0.5*(d1[k].lapse*d1[l].chi + d1[k].chi*d1[l].lapse - vars.h[k][l]*dlapse_dot_dchi);
        }
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM/2.0)*dlapse_dot_dchi;
    {
        FOR1(i)
        {
            tr_covd2lapse -= vars.chi*chrisvec[i]*d1[i].lapse;
            FOR1(j)
            {
                tr_covd2lapse += h_UU[i][j]*(vars.chi*d2[i][j].lapse + d1[i].lapse*d1[j].chi);
            }
        }
    }


    tensor<2, data_t> A_UU;
    {
        FOR2(i,j)
        {
            A_UU[i][j] = 0;
            FOR2(k,l)
            {
                A_UU[i][j] += h_UU[i][k]*h_UU[j][l]*vars.A[k][l];
            }
        }
    }

    data_t tr_ricci = 0;
    data_t tr_AA = 0;
    {
        FOR2(i,j)
        {
            tr_ricci += vars.chi*h_UU[i][j]*ricci[i][j];
            tr_AA += A_UU[i][j]*vars.A[i][j];
        }
    }

    {
        rhs.chi = advec.chi + (2.0/GR_SPACEDIM)*vars.chi*(vars.lapse*vars.K - divshift);
        FOR2(i,j)
        {
            rhs.h[i][j] = advec.h[i][j] - 2.0*vars.lapse*vars.A[i][j] - (2.0/GR_SPACEDIM)*vars.h[i][j]*divshift;
            FOR1(k)
            {
                rhs.h[i][j] += vars.h[k][i]*d1[j].shift[k] + vars.h[k][j]*d1[i].shift[k];
            }
        }
    }

    {
        tensor<2, data_t> Adot_TF_expr;
        FOR2(i,j)
        {
            Adot_TF_expr[i][j] = -covd2lapse[i][j] + vars.chi*vars.lapse*ricci[i][j];
        }

        data_t Adot_TF_trace = 0;
        FOR2(i,j)
        {
            Adot_TF_trace += h_UU[i][j]*Adot_TF_expr[i][j];
        }

        tensor<2, data_t> Adot_TF;
        FOR2(i,j)
        {
            Adot_TF[i][j] = Adot_TF_expr[i][j] - (1.0/GR_SPACEDIM)*Adot_TF_trace*vars.h[i][j];
        }

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
    }

#ifdef COVARIANTZ4
    data_t kappa1_lapse = m_params.kappa1;
#else
    data_t kappa1_lapse = m_params.kappa1 * vars.lapse;
#endif

    {
        data_t Thetadot = 0.5*vars.lapse*(tr_ricci - tr_AA + ((GR_SPACEDIM-1)/(double) GR_SPACEDIM)*vars.K*vars.K - 2*vars.Theta*vars.K) - 0.5*vars.Theta*kappa1_lapse*((GR_SPACEDIM+1) + m_params.kappa2*(GR_SPACEDIM-1)) - Z_dot_d1lapse;
        rhs.Theta = advec.Theta + Thetadot;

        rhs.K = advec.K + 2*Thetadot + vars.lapse*(tr_AA + (1.0/GR_SPACEDIM)*vars.K*vars.K) + kappa1_lapse*(1-m_params.kappa2)*vars.Theta + 2*Z_dot_d1lapse - tr_covd2lapse;
    }

    tensor<1, data_t> Gammadot;
    {
        FOR1(i)
        {
            Gammadot[i] = (2.0/GR_SPACEDIM)*(divshift*(chrisvec[i] + 2*m_params.kappa3*Z_over_chi[i]) - 2*vars.lapse*vars.K*Z_over_chi[i]) - 2*kappa1_lapse*Z_over_chi[i];
            FOR1(j)
            {
                Gammadot[i] += 2*h_UU[i][j]*(vars.lapse*d1[j].Theta - vars.Theta*d1[j].lapse)
                    - 2*A_UU[i][j]*d1[j].lapse
                    - vars.lapse*((2*(GR_SPACEDIM-1)/(double) GR_SPACEDIM)*h_UU[i][j]*d1[j].K + GR_SPACEDIM*A_UU[i][j]*d1[j].chi/chi_regularised)
                    - (chrisvec[j] + 2*m_params.kappa3*Z_over_chi[j])*d1[j].shift[i];

                FOR1(k)
                {
                    Gammadot[i] += 2*vars.lapse*chris[i][j][k]*A_UU[j][k]
                        + h_UU[j][k]*d2[j][k].shift[i]
                        + ((GR_SPACEDIM-2)/(double) GR_SPACEDIM)*h_UU[i][j]*d2[j][k].shift[k];
                }
            }
        }

        FOR1(i)
        {
            rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
        }
    }
    
    const data_t eta = 1;

    {
        rhs.lapse = m_params.lapse_advec_coeff*advec.lapse - 2*vars.lapse*(vars.K - 2*vars.Theta);
        FOR1(i)
        {
            rhs.shift[i] = m_params.shift_advec_coeff*advec.shift[i] + m_params.shift_gamma_coeff*vars.B[i];
            rhs.B[i] = m_params.shift_advec_coeff*advec.B[i] + (1 - m_params.shift_advec_coeff)*advec.Gamma[i] + Gammadot[i] - m_params.beta_driver*eta*vars.B[i];
        }
    }
}

template <class data_t>
void
CCZ4::local_vars(int idx, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = SIMDIFY<data_t>(m_in_ptr[i])[idx];
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::diff1(int idx, int stride, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far  = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        result[i] = (  weight_far  * in[idx - 2*stride]
                     - weight_near * in[idx - stride]
                     + weight_near * in[idx + stride]
                     - weight_far  * in[idx + 2*stride]) / m_dx;
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::diff2(int idx, int stride, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far   = 8.33333333333333333333e-2;
        data_t weight_near  = 1.33333333333333333333e+0;
        data_t weight_local = 2.50000000000000000000e+0;

        result[i] = (- weight_far   * in[idx - 2*stride]
                     + weight_near  * in[idx - stride]
                     - weight_local * in[idx]
                     + weight_near  * in[idx + stride]
                     - weight_far   * in[idx + 2*stride]) / (m_dx*m_dx);
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::mixed_diff2(int idx, int stride1, int stride2, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far_far   = 6.94444444444444444444e-3;
        data_t weight_near_far  = 5.55555555555555555556e-2;
        data_t weight_near_near = 4.44444444444444444444e-1;

        result[i] = (  weight_far_far   * in[idx - 2*stride1 - 2*stride2]
                     - weight_near_far  * in[idx - 2*stride1 - stride2]
                     + weight_near_far  * in[idx - 2*stride1 + stride2]
                     - weight_far_far   * in[idx - 2*stride1 + 2*stride2]

                     - weight_near_far  * in[idx - stride1 - 2*stride2]
                     + weight_near_near * in[idx - stride1 - stride2]
                     - weight_near_near * in[idx - stride1 + stride2]
                     + weight_near_far  * in[idx - stride1 + 2*stride2]

                     + weight_near_far  * in[idx + stride1 - 2*stride2]
                     - weight_near_near * in[idx + stride1 - stride2]
                     + weight_near_near * in[idx + stride1 + stride2]
                     - weight_near_far  * in[idx + stride1 + 2*stride2]

                     - weight_far_far   * in[idx + 2*stride1 - 2*stride2]
                     + weight_near_far  * in[idx + 2*stride1 - stride2]
                     - weight_near_far  * in[idx + 2*stride1 + stride2]
                     + weight_far_far   * in[idx + 2*stride1 + 2*stride2]) / (m_dx*m_dx);
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::advection(int idx, const tensor<1, data_t>& shift, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = 0;
    }

    for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
    {
        const int stride = m_stride[dim];
        const auto shift_val = shift[dim];
        const auto shift_positive = simd_compare_gt(shift_val, 0.0);

        for (int i = 0; i < c_NUM; ++i)
        {
            const auto in = SIMDIFY<data_t>(m_in_ptr[i]);
            const data_t in_left = in[idx - stride];
            const data_t in_centre = in[idx];
            const data_t in_right = in[idx + stride];
            
            data_t weight_0 = -2.50000000000000000000e-1;
            data_t weight_1 = -8.33333333333333333333e-1;
            data_t weight_2 = +1.50000000000000000000e+0;
            data_t weight_3 = -5.00000000000000000000e-1;
            data_t weight_4 = +8.33333333333333333333e-2;

            data_t upwind;
            upwind = shift_val * (  weight_0 * in_left
                                  + weight_1 * in_centre
                                  + weight_2 * in_right
                                  + weight_3 * in[idx + 2*stride]
                                  + weight_4 * in[idx + 3*stride]) / m_dx;

            data_t downwind;
            downwind = shift_val * (- weight_4 * in[idx - 3*stride]
                                    - weight_3 * in[idx - 2*stride]
                                    - weight_2 * in_left
                                    - weight_1 * in_centre
                                    - weight_0 * in_right) / m_dx;

            result[i] += simd_conditional(shift_positive, upwind , downwind);
        }
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::dissipation(int idx, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = 0;
    }

    for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
    {
        const int stride = m_stride[dim];
        for (int i = 0; i < c_NUM; ++i)
        {
            auto in = SIMDIFY<data_t>(m_in_ptr[i]);

            data_t weight_vfar  = 1.56250e-2;
            data_t weight_far   = 9.37500e-2;
            data_t weight_near  = 2.34375e-1;
            data_t weight_local = 3.12500e-1;

            result[i] += ( weight_vfar   * in[idx - 3*stride]
                          - weight_far   * in[idx - 2*stride]
                          + weight_near  * in[idx - stride]
                          - weight_local * in[idx]
                          + weight_near  * in[idx + stride]
                          - weight_far   * in[idx + 2*stride]
                          + weight_vfar  * in[idx + 3*stride]) / m_dx;
        }
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out)
{
    out.chi      = in[c_chi];
    out.h[0][0]  = in[c_h11];
    out.h[0][1]  = in[c_h12];
    out.h[0][2]  = in[c_h13];
    out.h[1][1]  = in[c_h22];
    out.h[1][2]  = in[c_h23];
    out.h[2][2]  = in[c_h33];

    out.h[1][0] = out.h[0][1];
    out.h[2][0] = out.h[0][2];
    out.h[2][1] = out.h[1][2];

    out.K        = in[c_K];
    out.A[0][0]  = in[c_A11];
    out.A[0][1]  = in[c_A12];
    out.A[0][2]  = in[c_A13];
    out.A[1][1]  = in[c_A22];
    out.A[1][2]  = in[c_A23];
    out.A[2][2]  = in[c_A33];

    out.A[1][0] = out.A[0][1];
    out.A[2][0] = out.A[0][2];
    out.A[2][1] = out.A[1][2];

    out.Gamma[0] = in[c_Gamma1];
    out.Gamma[1] = in[c_Gamma2];
    out.Gamma[2] = in[c_Gamma3];

    out.Theta    = in[c_Theta];

    out.lapse    = in[c_lapse];
    out.shift[0] = in[c_shift1];
    out.shift[1] = in[c_shift2];
    out.shift[2] = in[c_shift3];

    out.B[0]     = in[c_B1];
    out.B[1]     = in[c_B2];
    out.B[2]     = in[c_B3];
}
