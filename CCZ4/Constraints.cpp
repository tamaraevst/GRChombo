#include "Constraints.hpp"


template <class data_t>
void
Constraints::compute(int x, int y, int z)
{
    const int idx = m_driver.m_stride[2]*(z-m_driver.m_in_lo[2]) + m_driver.m_stride[1]*(y-m_driver.m_in_lo[1]) + (x-m_driver.m_in_lo[0]);

    vars_t<data_t> vars;
    local_vars(idx, vars);

    vars_t<data_t> d1[CH_SPACEDIM];
    for (int i = 0; i < 3; ++i)
    {
        diff1(idx, m_driver.m_stride[i], d1[i]);
    }

    vars_t<data_t> d2[CH_SPACEDIM][CH_SPACEDIM];

    // Repeated derivatives
    for (int i = 0; i < 3; ++i)
    {
        diff2(idx, m_driver.m_stride[i], d2[i][i]);
    }

    // Mixed derivatives
    mixed_diff2(idx, m_driver.m_stride[1], m_driver.m_stride[0], d2[0][1]);
    mixed_diff2(idx, m_driver.m_stride[2], m_driver.m_stride[0], d2[0][2]);
    mixed_diff2(idx, m_driver.m_stride[2], m_driver.m_stride[1], d2[1][2]);

    d2[1][0] = d2[0][1];
    d2[2][0] = d2[0][2];
    d2[2][1] = d2[1][2];

    constraint_equations(vars, d1, d2);
}

template <class data_t>
void
Constraints::constraint_equations(vars_t<data_t> &vars,
      const vars_t<data_t> (&d1)[CH_SPACEDIM],
      const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]
      )
{
   const data_t chi_regularised = simd_max(1e-6, vars.chi);

   auto h_UU = compute_inverse_metric(vars);
   chris_t<data_t> chris(vars, d1, h_UU);

   ricci_t<data_t> ricci(vars, d2, d1, chris, h_UU);

   auto A_UU       = raise(vars.A, h_UU);
   data_t tr_AA    = trace(vars.A, A_UU);

   vars.HamC = ricci.scalar + (GR_SPACEDIM-2.)*pow(vars.K,2)/(GR_SPACEDIM-1.) - tr_AA;

   tensor<1,data_t> momC;

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
      momC[i] = - (GR_SPACEDIM-2.)*d1[i].K/(GR_SPACEDIM-1.);
   }
   FOR3(i,j,k)
   {
      momC[i] += h_UU(j,k)*(covd_A[k][j][i] - (GR_SPACEDIM-1.) * vars.A[i][j] * d1[k] / (2 * chi_regularised));
   }
   vars.MomC1 = momC[0]; //Write to grid functions
   vars.MomC2 = momC[1];
   vars.MomC3 = momC[2];
}
