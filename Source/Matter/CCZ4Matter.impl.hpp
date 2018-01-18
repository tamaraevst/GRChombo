// Last edited K Clough 16.02.17

#if !defined(CCZ4MATTER_HPP_)
#error "This file should only be included through CCZ4Matter.hpp"
#endif

#ifndef CCZ4MATTER_IMPL_HPP_
#define CCZ4MATTER_IMPL_HPP_

#define COVARIANTZ4
#include "MiscUtils.hpp"

template <class matter_t>
CCZ4Matter<matter_t>::CCZ4Matter(matter_t a_matter, params_t params, double dx,
                                 double sigma, int formulation, double G_Newton)
    : CCZ4(params, dx, sigma, formulation, 0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void CCZ4Matter<matter_t>::compute(Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto matter_vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

    // Call CCZ4 RHS - work out RHS without matter, no dissipation
    Vars<data_t> matter_rhs;
    rhs_equation(matter_rhs, matter_vars, d1, d2, advec);

    // add RHS matter terms from EM tensor
    add_EMTensor_rhs(matter_rhs, matter_vars, d1);

    // add evolution of matter fields themselves
    my_matter.add_matter_rhs(matter_rhs, matter_vars, d1, d2, advec);

    // Add dissipation to all terms
    m_deriv.add_dissipation(matter_rhs, current_cell, m_sigma);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(matter_rhs);
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t>
template <class data_t>
void CCZ4Matter<matter_t>::add_EMTensor_rhs(
    Vars<data_t> &matter_rhs, const Vars<data_t> &matter_vars,
    const Vars<tensor<1, data_t>> &d1) const
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(matter_vars.h);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Calculate elements of the decomposed stress energy tensor
    const auto emtensor =
        my_matter.compute_emtensor(matter_vars, d1, h_UU, chris.ULL);

    // Update RHS for K and Theta depending on formulation
    if (m_formulation == USE_BSSN)
    {
        matter_rhs.K += 4.0 * M_PI * m_G_Newton * matter_vars.lapse *
                        (emtensor.S + emtensor.rho);
        matter_rhs.Theta += 0.0;
    }
    else
    {
        matter_rhs.K += 4.0 * M_PI * m_G_Newton * matter_vars.lapse *
                        (emtensor.S - 3 * emtensor.rho);
        matter_rhs.Theta +=
            -8.0 * M_PI * m_G_Newton * matter_vars.lapse * emtensor.rho;
    }

    // Update RHS for other variables
    tensor<2, data_t> Sij_TF = emtensor.Sij;
    make_trace_free(Sij_TF, matter_vars.h, h_UU);

    FOR2(i, j)
    {
        matter_rhs.A[i][j] += -8.0 * M_PI * m_G_Newton * matter_vars.chi *
                              matter_vars.lapse * Sij_TF[i][j];
    }

    FOR1(i)
    {
        data_t matter_term_Gamma = 0.0;
        FOR1(j)
        {
            matter_term_Gamma += -16.0 * M_PI * m_G_Newton * matter_vars.lapse *
                                 h_UU[i][j] * emtensor.Si[j];
        }

        matter_rhs.Gamma[i] += matter_term_Gamma;
        matter_rhs.B[i] += matter_term_Gamma;
    }
}

#endif /* CCZ4MATTER_IMPL_HPP_ */
