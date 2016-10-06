#if !defined(RELAXATIONCHI_HPP_)
#error "This file should only be included through RelaxationChi.hpp"
#endif

#ifndef RELAXATIONCHI_IMPL_HPP_
#define RELAXATIONCHI_IMPL_HPP_


RelaxationChi::RelaxationChi(const FABDriverBase& driver, double dx, double relaxspeed) :
    m_relaxspeed (relaxspeed),
    m_driver (driver),
    m_deriv (dx, m_driver)
{}

template <class data_t>
void
RelaxationChi::compute(int ix, int iy, int iz)
{
    vars_t<data_t> vars;
    m_driver.local_vars(vars);  //copy data from chombo gridpoint into local variables

    vars_t< tensor<1, data_t> > d1;
    FOR1(idir) m_deriv.diff1(d1, idir); //work out first derivatives of variables on grid

    vars_t< tensor<2,data_t> > d2;
    // Repeated derivatives
    FOR1(idir) m_deriv.diff2(d2, idir);  //work out second derivatives of variables on grid
    // Mixed derivatives
    // Note: no need to symmetrise explicitely, this is done in mixed_diff2
    m_deriv.mixed_diff2(d2, 1, 0);
    m_deriv.mixed_diff2(d2, 2, 0);
    m_deriv.mixed_diff2(d2, 2, 1);

    vars_t<data_t> advec;
    advec.assign(0.);
    FOR1(idir) m_deriv.add_advection(advec, vars.shift[idir], idir);

    vars_t<data_t> rhs = rhs_equation(vars, d1, d2, advec); //work out RHS including advection

//    No dissipation in relaxation
//    FOR1(idir) m_deriv.add_dissipation(rhs, m_sigma, idir); //add dissipation to RHS

    //Write the rhs into the output FArrayBox
    m_driver.store_vars(rhs);
}

template <class data_t>
auto
RelaxationChi::rhs_equation(const vars_t<data_t> &vars,
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

    //Calculate elements of the decomposed stress energy tensor and ricci tensor
    auto emtensor =  CCZ4EMTensorSF::compute_emtensor_SF(vars, d1, h_UU, chris.ULL, advec.phi);
    auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

    auto A_UU       = TensorAlgebra::raise_all(vars.A, h_UU);
    data_t tr_AA    = TensorAlgebra::compute_trace(vars.A, A_UU);

    rhs.chi = - m_relaxspeed*(ricci.scalar + (GR_SPACEDIM-2.)*vars.K*vars.K/(GR_SPACEDIM-1.) - tr_AA - 16.0*M_PI*emtensor.rho);

    return rhs;
}

template <class data_t>
RelaxationChi::vars_t<data_t>::vars_t()
{
    //Define the mapping from components of chombo grid to elements in vars_t.
    //This allows to read/write data from the chombo grid into local
    //variables in vars_t (which only exist for the current cell).
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

    define_enum_mapping(c_Theta, Theta);

    define_enum_mapping(c_lapse, lapse);
    define_enum_mapping(c_shift1, shift[0]);
    define_enum_mapping(c_shift2, shift[1]);
    define_enum_mapping(c_shift3, shift[2]);

    define_enum_mapping(c_B1, B[0]);
    define_enum_mapping(c_B2, B[1]);
    define_enum_mapping(c_B3, B[2]);

    define_enum_mapping(c_phi, phi);
    define_enum_mapping(c_PiM, PiM);


}

#endif /* RELAXATIONCHI_IMPL_HPP_ */
