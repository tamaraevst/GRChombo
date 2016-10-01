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
    vars_t<data_t> vars;
    m_driver.local_vars(vars);

    vars_t< tensor<1, data_t> > d1;
    FOR1(idir) m_deriv.diff1(d1, idir);

    vars_t< tensor<2,data_t> > d2;
    // Repeated derivatives
    FOR1(idir) m_deriv.diff2(d2, idir);
    // Mixed derivatives
    // Note: no need to symmetrise explicitely, this is done in mixed_diff2
    m_deriv.mixed_diff2(d2, 1, 0);
    m_deriv.mixed_diff2(d2, 2, 0);
    m_deriv.mixed_diff2(d2, 2, 1);

    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    //Write the rhs into the output FArrayBox
    m_driver.store_vars(out.Ham, c_Ham);
    m_driver.store_vars(out.Mom[0], c_Mom1);
    m_driver.store_vars(out.Mom[1], c_Mom2);
    m_driver.store_vars(out.Mom[2], c_Mom3);
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
      FOR1(l)
      {
          covd_A[i][j][k] += - chris.ULL[l][i][j]*vars.A[l][k] - chris.ULL[l][i][k]*vars.A[l][j];
      }
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
    define_enum_mapping(c_chi, chi);

    define_enum_mapping(c_h11, h[0][0]);
    define_enum_mapping(c_h12, h[0][1]);
    define_enum_mapping(c_h12, h[1][0]);
    define_enum_mapping(c_h13, h[0][2]);
    define_enum_mapping(c_h13, h[2][0]);
    define_enum_mapping(c_h22, h[1][1]);
    define_enum_mapping(c_h23, h[1][2]);
    define_enum_mapping(c_h23, h[2][1]);
    define_enum_mapping(c_h33, h[2][2]);

    define_enum_mapping(c_K, K);

    define_enum_mapping(c_A11, A[0][0]);
    define_enum_mapping(c_A12, A[0][1]);
    define_enum_mapping(c_A12, A[1][0]);
    define_enum_mapping(c_A13, A[0][2]);
    define_enum_mapping(c_A13, A[2][0]);
    define_enum_mapping(c_A22, A[1][1]);
    define_enum_mapping(c_A23, A[1][2]);
    define_enum_mapping(c_A23, A[2][1]);
    define_enum_mapping(c_A33, A[2][2]);

    define_enum_mapping(c_Gamma1, Gamma[0]);
    define_enum_mapping(c_Gamma2, Gamma[1]);
    define_enum_mapping(c_Gamma3, Gamma[2]);
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
