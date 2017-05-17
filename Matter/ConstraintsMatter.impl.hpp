// Last edited K Clough 16.02.17

#if !defined(CONSTRAINTSMATTER_HPP_)
#error "This file should only be included through ConstraintsMatter.hpp"
#endif

#ifndef CONSTRAINTSMATTER_IMPL_HPP_
#define CONSTRAINTSMATTER_IMPL_HPP_

template <class matter_t>
ConstraintsMatter<matter_t>::ConstraintsMatter(
    const FABDriverBase& driver,
    const matter_t a_matter,
    double dx,
    double G_Newton)
    : Constraints(driver, dx, 0.0 /*No cosmological constant*/),
      my_matter (a_matter),
      m_G_Newton (G_Newton) {}


template <class matter_t>
template <class data_t>
void ConstraintsMatter<matter_t>::compute(Cell current_cell)
{
    //Calculate non matter contributions to Constraints
    Vars<data_t> vars;
    m_driver.local_vars(vars, current_cell);

    //Calculate first derivatives
    Vars< tensor<1, data_t> > d1;
    FOR1(idir) m_deriv.diff1(d1, current_cell, idir);

    //Calculate second derivatives
    Vars< tensor<2,data_t> > d2;
    // Repeated derivatives
    FOR1(idir) m_deriv.diff2(d2, current_cell, idir);
    // Mixed derivatives
    // Note: no need to symmetrise explicitely, this is done in mixed_diff2
    m_deriv.mixed_diff2(d2, current_cell, 1, 0);
    m_deriv.mixed_diff2(d2, current_cell, 2, 0);
    m_deriv.mixed_diff2(d2, current_cell, 2, 1);

    // Get the non matter terms
    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    //Calculate EM Tensor and add matter terms, need advection and geometric objects
    Vars<data_t> advec;
    advec.assign(0.);
    FOR1(idir) m_deriv.add_advection(advec, current_cell, vars.shift[idir], idir);

    // Inverse metric and Christoffel symbol
    auto h_UU = TensorAlgebra::compute_inverse(vars.h);
    auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Energy Momentum tensor
    auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL, advec);

    //Hamiltonain constraint
    out.Ham += -16.0*M_PI*m_G_Newton*emtensor.rho;

    //Momentum constraints
    FOR1(i)
    {
        out.Mom[i] += -8.0*M_PI*m_G_Newton*emtensor.Si[i];
    }

    //Write the rhs into the output FArrayBox
    m_driver.store_vars(out.Ham, current_cell, c_Ham);
    m_driver.store_vars(out.Mom[0], current_cell, c_Mom1);
    m_driver.store_vars(out.Mom[1], current_cell, c_Mom2);
    m_driver.store_vars(out.Mom[2], current_cell, c_Mom3);
}

#endif /* CONSTRAINTSMATTER_IMPL_HPP_ */
