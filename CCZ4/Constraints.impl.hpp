#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

Constraints::Constraints(const FABDriverBase& driver, double dx) :
    m_driver (driver),
    m_deriv (dx, m_driver)
{}

template <class data_t>
void
Constraints::compute(int ix, int iy, int iz)
{
    idx_t<data_t> idx = m_driver.in_idx(ix, iy, iz);

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

    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    idx_t<data_t> out_idx = m_driver.out_idx(ix, iy, iz);
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Ham])[out_idx]  = out.Ham;
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom1])[out_idx] = out.Mom[0];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom2])[out_idx] = out.Mom[1];
    SIMDIFY<data_t>(m_driver.m_out_ptr[c_Mom3])[out_idx] = out.Mom[2];
}

template <class data_t>
auto
Constraints::constraint_equations(
      vars_t<data_t> &vars,
      const vars_t< tensor<1,data_t> >& d1,
      const vars_t< tensor<2,data_t> >& d2
) -> constraints_t<data_t>
{
   constraints_t<data_t> out;

   const data_t chi_regularised = simd_max(1e-6, vars.chi);

   auto h_UU = TensorAlgebra::compute_inverse(vars.h);
   auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

   auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

   auto A_UU       = TensorAlgebra::raise_all(vars.A, h_UU);
   data_t tr_AA    = TensorAlgebra::compute_trace(vars.A, A_UU);

   out.Ham = ricci.scalar + (GR_SPACEDIM-2.)*vars.K*vars.K/(GR_SPACEDIM-1.) - tr_AA;

   tensor<2,data_t> covd_A[CH_SPACEDIM];
   FOR3(i,j,k)
   {
      covd_A[i][j][k] = d1.A[j][k][i];
   }
   FOR4(i,j,k,l)
   {
      covd_A[i][j][k] += - chris.ULL[l][i][j]*vars.A[l][k] - chris.ULL[l][i][k]*vars.A[l][j];
   }

   FOR1(i)
   {
      out.Mom[i] = - (GR_SPACEDIM-1.)*d1.K[i]/GR_SPACEDIM;
   }
   FOR3(i,j,k)
   {
      out.Mom[i] += h_UU[j][k]*(covd_A[k][j][i] - GR_SPACEDIM * vars.A[i][j] * d1.chi[k] / (2 * chi_regularised));
   }

   return out;
}

template <class data_t>
Constraints::vars_t<data_t>::vars_t()
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
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
