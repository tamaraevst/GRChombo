// Last edited K Clough 15.02.17

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template<typename> class vars_t>
auto ScalarField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars,
    const vars_t< tensor<1,data_t> >& d1,
    const tensor<2, data_t> &h_UU,
    const tensor<3, data_t> &chris_ULL) -> emtensor_t<data_t> {

    emtensor_t<data_t> out;

    //set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    //compute potential
    my_potential.compute_potential(V_of_phi, dVdphi, vars.phi);

    // Useful quantity Vt
    data_t Vt = - vars.Pi * vars.Pi + 2.0*V_of_phi;
    FOR2(i,j)
    {
        Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i,j)
    {
        out.Sij[i][j] = - 0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij , h_UU);

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] = - d1.phi[i]*vars.Pi;
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi*vars.Pi + 0.5*Vt;

    return out;
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template<typename> class vars_t, template<typename> class diff2_vars_t>
void ScalarField<potential_t>::add_matter_rhs(
    vars_t<data_t> &total_rhs,
    const vars_t<data_t> &vars,
    const vars_t< tensor<1,data_t> >& d1,
    const diff2_vars_t< tensor<2,data_t> >& d2,
    const vars_t<data_t> &advec) {

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    //set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;

    //compute potential
    my_potential.compute_potential(V_of_phi, dVdphi, vars.phi);

    //evolution equations for scalar field and (minus) its conjugate momentum
    total_rhs.phi = vars.lapse*vars.Pi + advec.phi;
    total_rhs.Pi = vars.lapse*(vars.K*vars.Pi - dVdphi) + advec.Pi;

    FOR2(i,j)
    {
        //includes non conformal parts of chris not included in chris_ULL
        total_rhs.Pi += h_UU[i][j]*( - 0.5*d1.chi[j]*vars.lapse*d1.phi[i]
                                       + vars.chi*vars.lapse*d2.phi[i][j]
                                       + vars.chi*d1.lapse[i]*d1.phi[j]  );
        FOR1(k)
        {
            total_rhs.Pi += - vars.chi*vars.lapse*h_UU[i][j]*chris.ULL[k][i][j]*d1.phi[k];
        }
    }
}

#endif /* SCALARFIELD_IMPL_HPP_ */
