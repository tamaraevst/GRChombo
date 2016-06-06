#include "Constraints.hpp"

Constraints::Constraints(double dx, const FABDriverBase& driver) :
    m_dx (dx),
    m_driver (driver)
{}

template <class data_t>
void
Constraints::compute(int x, int y, int z)
{
    const int idx = m_driver.m_stride[2]*(z-m_driver.m_in_lo[2]) + m_driver.m_stride[1]*(y-m_driver.m_in_lo[1]) + (x-m_driver.m_in_lo[0]);

    vars_t<data_t> vars;
    {
       data_t varsArr[c_NUM];
       m_driver.local_vars(idx, varsArr);
       demarshall(varsArr, vars);
    }

    vars_t<data_t> d1[CH_SPACEDIM];
    for (int i = 0; i < 3; ++i)
    {
       {
          data_t varsArr[c_NUM];
          m_driver.diff1(idx, m_driver.m_stride[i], m_dx, varsArr);
          demarshall(varsArr, d1[i]);
       }
    }

    vars_t<data_t> d2[CH_SPACEDIM][CH_SPACEDIM];

    // Repeated derivatives
    for (int i = 0; i < 3; ++i)
    {
       {
          data_t varsArr[c_NUM];
          m_driver.diff2(idx, m_driver.m_stride[i], m_dx, varsArr);
          demarshall(varsArr, d2[i][i]);
       }
    }

    // Mixed derivatives
    {
       data_t varsArr[c_NUM];
       m_driver.mixed_diff2(idx, m_driver.m_stride[1], m_driver.m_stride[0], m_dx, varsArr);
       demarshall(varsArr, d2[0][1]);
       m_driver.mixed_diff2(idx, m_driver.m_stride[2], m_driver.m_stride[0], m_dx, varsArr);
       demarshall(varsArr, d2[0][2]);
       m_driver.mixed_diff2(idx, m_driver.m_stride[2], m_driver.m_stride[1], m_dx, varsArr);
       demarshall(varsArr, d2[1][2]);
    }

    d2[1][0] = d2[0][1];
    d2[2][0] = d2[0][2];
    d2[2][1] = d2[1][2];

    constraint_equations(vars, d1, d2);

    // TODO: I really do not like this, but cannot think of a better way to do it yet...
    const int out_idx = m_driver.m_out_stride[2]*(z-m_driver.m_out_lo[2]) + m_driver.m_out_stride[1]*(y-m_driver.m_out_lo[1]) + (x-m_driver.m_out_lo[0]);
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Ham])[out_idx]      = vars.Ham;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom1])[out_idx]     = vars.Mom[0];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom2])[out_idx]     = vars.Mom[1];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom3])[out_idx]     = vars.Mom[2];
}

template <class data_t>
void
Constraints::constraint_equations(vars_t<data_t> &vars,
      const vars_t<data_t> (&d1)[CH_SPACEDIM],
      const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]
      )
{
   const data_t chi_regularised = simd_max(1e-6, vars.chi);

   auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);
   auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

   auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

   auto A_UU       = raiseAll(vars.A, h_UU);
   data_t tr_AA    = compute_trace(vars.A, A_UU);

   vars.Ham = ricci.scalar + (GR_SPACEDIM-2.)*vars.K*vars.K/(GR_SPACEDIM-1.) - tr_AA;

   tensor<2,data_t> covd_A[CH_SPACEDIM];
   FOR3(i,j,k)
   {
      covd_A[i][j][k] = d1[i].A[j][k];
   }
   FOR4(i,j,k,l)
   {
      covd_A[i][j][k] += - chris.ULL[l][i][j]*vars.A[l][k] - chris.ULL[l][i][k]*vars.A[l][j];
   }

   FOR1(i)
   {
      vars.Mom[i] = - (GR_SPACEDIM-1.)*d1[i].K/GR_SPACEDIM;
   }
   FOR3(i,j,k)
   {
      vars.Mom[i] += h_UU[j][k]*(covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] * d1[k].chi / (2 * chi_regularised));
   }
}

template <class data_t>
void
Constraints::demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out)
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

    out.Ham    = in[c_Ham];

    out.Mom[0]     = in[c_Mom1];
    out.Mom[1]     = in[c_Mom2];
    out.Mom[2]     = in[c_Mom3];
}
