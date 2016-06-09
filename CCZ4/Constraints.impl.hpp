#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

Constraints::Constraints(double dx, const FABDriverBase& driver) :
    m_driver (driver),
    m_deriv (dx, m_driver)
{}

template <class data_t>
void
Constraints::compute(int x, int y, int z)
{
    idx_t<data_t> idx = m_driver.in_idx(x, y, z);

    vars_t<data_t> vars = m_driver.local_vars(idx);

    vars_t<data_t> d1[CH_SPACEDIM];
    for (int dir = 0; dir < 3; ++dir)
    {
      d1[dir] = m_deriv.diff1(idx, dir), d1[dir];
    }

    vars_t<data_t> d2[CH_SPACEDIM][CH_SPACEDIM];

    // Repeated derivatives
    for (int dir = 0; dir < 3; ++dir)
    {
      d2[dir][dir] = m_deriv.diff2(idx, dir), d2[dir][dir];
    }

    // Mixed derivatives
    d2[0][1] = m_deriv.mixed_diff2<data_t>(idx, 1, 0);
    d2[0][2] = m_deriv.mixed_diff2<data_t>(idx, 2, 0);
    d2[1][2] = m_deriv.mixed_diff2<data_t>(idx, 2, 1);
    
    d2[1][0] = d2[0][1];
    d2[2][0] = d2[0][2];
    d2[2][1] = d2[1][2];

    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Ham])[out_idx]  = out.Ham;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom1])[out_idx] = out.Mom[0];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom2])[out_idx] = out.Mom[1];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom3])[out_idx] = out.Mom[2];
}

template <class data_t>
auto
Constraints::constraint_equations(
      vars_t<data_t> &vars,
      const vars_t<data_t> (&d1)[CH_SPACEDIM],
      const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]
) -> constraints_t<data_t>
{
   constraints_t<data_t> out;

   const data_t chi_regularised = simd_max(1e-6, vars.chi);

   auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);
   auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

   auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

   auto A_UU       = raise_all(vars.A, h_UU);
   data_t tr_AA    = compute_trace(vars.A, A_UU);

   out.Ham = ricci.scalar + (GR_SPACEDIM-2.)*vars.K*vars.K/(GR_SPACEDIM-1.) - tr_AA;

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
      out.Mom[i] = - (GR_SPACEDIM-1.)*d1[i].K/GR_SPACEDIM;
   }
   FOR3(i,j,k)
   {
      out.Mom[i] += h_UU[j][k]*(covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] * d1[k].chi / (2 * chi_regularised));
   }

   return out;
}

template <class data_t> template <class arr_t>
Constraints::vars_t<data_t>::vars_t(const arr_t& in)
{
    chi      = in[c_chi];
    h[0][0]  = in[c_h11];
    h[0][1]  = in[c_h12];
    h[0][2]  = in[c_h13];
    h[1][1]  = in[c_h22];
    h[1][2]  = in[c_h23];
    h[2][2]  = in[c_h33];

    h[1][0] = h[0][1];
    h[2][0] = h[0][2];
    h[2][1] = h[1][2];

    K        = in[c_K];
    A[0][0]  = in[c_A11];
    A[0][1]  = in[c_A12];
    A[0][2]  = in[c_A13];
    A[1][1]  = in[c_A22];
    A[1][2]  = in[c_A23];
    A[2][2]  = in[c_A33];

    A[1][0] = A[0][1];
    A[2][0] = A[0][2];
    A[2][1] = A[1][2];

    Gamma[0] = in[c_Gamma1];
    Gamma[1] = in[c_Gamma2];
    Gamma[2] = in[c_Gamma3];
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
