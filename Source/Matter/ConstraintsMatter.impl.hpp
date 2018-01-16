// Last edited K Clough 16.02.17

#if !defined(CONSTRAINTSMATTER_HPP_)
#error "This file should only be included through ConstraintsMatter.hpp"
#endif

#ifndef CONSTRAINTSMATTER_IMPL_HPP_
#define CONSTRAINTSMATTER_IMPL_HPP_

template <class matter_t>
ConstraintsMatter<matter_t>::ConstraintsMatter(
    const matter_t a_matter,
    double dx,
    double G_Newton)
    : Constraints(dx, 0.0 /*No cosmological constant*/),
      my_matter (a_matter),
      m_G_Newton (G_Newton) {}

template <class matter_t>
template <class data_t>
void ConstraintsMatter<matter_t>::compute(Cell<data_t> current_cell)
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    // Get the non matter terms for the constraints
    constraints_t<data_t> out = constraint_equations(vars, d1, d2);

    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Energy Momentum tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    //Hamiltonain constraint
    out.Ham += -16.0*M_PI*m_G_Newton*emtensor.rho;

    //Momentum constraints
    FOR1(i)
    {
        out.Mom[i] += -8.0*M_PI*m_G_Newton*emtensor.Si[i];
    }

    //Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Ham, c_Ham);
    current_cell.store_vars(out.Mom, GRInterval<c_Mom1,c_Mom3>());
}

#endif /* CONSTRAINTSMATTER_IMPL_HPP_ */
