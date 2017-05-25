// Last edited K Clough 16.02.17

#if !defined(RELAXATIONCHI_HPP_)
#error "This file should only be included through RelaxationChi.hpp"
#endif

#ifndef RELAXATIONCHI_IMPL_HPP_
#define RELAXATIONCHI_IMPL_HPP_

template <class matter_t>
RelaxationChi<matter_t>::RelaxationChi(
    matter_t a_matter,
    double dx,
    double relax_speed,
    double G_Newton)
    : my_matter (a_matter),
      m_relax_speed (relax_speed),
      m_G_Newton (G_Newton),
      m_deriv (dx) {}

template <class matter_t>
template <class data_t>
void RelaxationChi<matter_t>::compute(Cell current_cell) {

    //copy data from chombo gridpoint into local variables
    Vars<data_t> vars;
    current_cell.local_vars(vars);

    //work out first derivatives of variables on grid
    Vars< tensor<1, data_t> > d1;
    FOR1(idir) m_deriv.diff1(d1, current_cell, idir);

    //work out second derivatives of variables on grid
    Vars< tensor<2,data_t> > d2;
    // Repeated derivatives
    FOR1(idir) m_deriv.diff2(d2, current_cell, idir);
    // Mixed derivatives
    // Note: no need to symmetrise explicitely, this is done in mixed_diff2
    m_deriv.mixed_diff2(d2, current_cell, 1, 0);
    m_deriv.mixed_diff2(d2, current_cell, 2, 0);
    m_deriv.mixed_diff2(d2, current_cell, 2, 1);

    // Calculate advection components
    Vars<data_t> advec;
    advec.assign(0.);
    FOR1(idir) m_deriv.add_advection(advec, current_cell, vars.shift[idir], idir);

    //work out RHS including advection
    Vars<data_t> rhs;
    rhs.assign(0.); //All components that are not explicitly set in rhs_equation are 0
    rhs_equation(rhs, vars, d1, d2, advec);

    //    No dissipation in relaxation for now but may add it
    //    FOR1(idir) m_deriv.add_dissipation(rhs, m_sigma, idir);

    //Write the rhs into the output FArrayBox
    current_cell.store_vars(rhs);
}

template <class matter_t>
template <class data_t>
void RelaxationChi<matter_t>::rhs_equation(
    Vars<data_t> &rhs,
    const Vars<data_t> &vars,
    const Vars< tensor<1,data_t> >& d1,
    const Vars< tensor<2,data_t> >& d2,
    const Vars<data_t> &advec) {

    using namespace TensorAlgebra;

    auto h_UU = compute_inverse_sym(vars.h);
    auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    //Calculate elements of the decomposed stress energy tensor and ricci tensor
    auto emtensor =  my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL, advec);
    auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);
    auto A_UU       = raise_all(vars.A, h_UU);
    data_t tr_AA    = compute_trace(vars.A, A_UU);

    //Calculate the relaxation RHS for chi, all other vars RHS zero
    //Could have called ConstraintsMatter here, but it is hardly worth it
    rhs.chi =  m_relax_speed*(ricci.scalar+(GR_SPACEDIM-1.)*vars.K*vars.K/GR_SPACEDIM
                           - tr_AA - 16.0*M_PI*m_G_Newton*emtensor.rho);
}

#endif /* RELAXATIONCHI_IMPL_HPP_ */
