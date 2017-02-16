// Last edited K Clough 15.02.17

#if !defined(SCALARFIELD_HPP_)
#error "This file should only be included through ScalarField.hpp"
#endif

#ifndef SCALARFIELD_IMPL_HPP_
#define SCALARFIELD_IMPL_HPP_

#define COVARIANTZ4

// Calculate the stress energy tensor elements
template <class data_t>
auto ScalarField::compute_emtensor(
    const Vars<data_t> &vars,
    const Vars< tensor<1,data_t> >& d1,
    const tensor<2, data_t> &h_UU,
    const tensor<3, data_t> &chris_ULL,
    const Vars<data_t> &advec) -> emtensor_t<data_t> {

  emtensor_t<data_t> out;

  // Find the potential and its gradient in terms of phi
  potential_t<data_t> potential = compute_potential(vars.phi);

  // Some useful quantities
  data_t Vt = - vars.Pi * vars.Pi + 2.0*potential.V_of_phi;
  FOR2(i,j)
  {
    Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
  }

  data_t dphidt2;
  dphidt2 = (vars.lapse*vars.Pi + advec.phi)*(vars.lapse*vars.Pi + advec.phi);

  tensor<1, data_t> T_i; // The T(i,3) components of the 4d stress energy tensor
  FOR1(i)
  {
    T_i[i] = (d1.phi[i] * (vars.Pi*vars.lapse + advec.phi));

    FOR1(j)
    {
      T_i[i] += -0.5*Vt*vars.h[i][j]*vars.shift[j]/vars.chi;
    }
  }

  // Calculate components of EM Tensor
  // S_ij = T_ij
  FOR2(i,j)
  {
    out.Sij[i][j] = -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
  }

  // S = Tr_S_ij
  out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij , h_UU);

  // S_i (note lower index) = n^a T_a0
  FOR1(i)
  {
    out.Si[i] = - T_i[i]/vars.lapse;

    FOR1(j)
    {
      out.Si[i] += vars.shift[j]/vars.lapse * out.Sij[i][j];
    }
  }

  auto lapse_squared = vars.lapse*vars.lapse;
  // rho = n^a n^b T_ab
  out.rho = dphidt2/lapse_squared + 0.5*Vt;
  FOR2(i,j)
  {
    out.rho += (-0.5*Vt*vars.h[i][j]/vars.chi + out.Sij[i][j])
              * vars.shift[i]*vars.shift[j]/lapse_squared;
  }
  FOR1(i)
  {
    out.rho += - 2.0*vars.shift[i]*T_i[i]/lapse_squared;
  }

  return out;
}

/// Set the potential function for the scalar field here
template <class data_t>
auto ScalarField::compute_potential(const data_t phi_here) -> potential_t<data_t> {

  potential_t<data_t> out;

  //The potential value at phi
  out.V_of_phi = m_params.scalar_mass*phi_here; // e.g. m^2 phi^2

  //The potential gradient at phi
  out.dVdphi = m_params.scalar_mass;  //  e.g. 2 m^2 phi

  return out;
}

// Sums all contributions the the RHS, including matter terms
template <class data_t>
void ScalarField::add_matter_rhs(
    Vars<data_t> &total_rhs,
    const Vars<data_t> &vars,
    const Vars< tensor<1,data_t> >& d1,
    const Vars< tensor<2,data_t> >& d2,
    const Vars<data_t> &advec) {

  using namespace TensorAlgebra;

  auto h_UU = compute_inverse(vars.h);
  auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);
  potential_t<data_t> potential = compute_potential(vars.phi);

  //evolution equations for scalar field and (minus) its conjugate momentum
  total_rhs.phi = vars.lapse*vars.Pi + advec.phi;
  total_rhs.Pi = vars.lapse*(vars.K*vars.Pi - potential.dVdphi) + advec.Pi;

  FOR2(i,j)
  {
    //includes non conformal parts of chris not included in chris_ULL
    total_rhs.Pi += h_UU[i][j]*( - 0.5*d1.chi[j]*vars.lapse*d1.phi[i]
                                  + vars.chi*vars.lapse*d2.phi[i][j]
                                  + vars.chi*d1.lapse[i]*d1.phi[j]   );
    FOR1(k)
    {
      total_rhs.Pi += - vars.chi * vars.lapse * h_UU[i][j]
                          * chris.ULL[k][i][j] * d1.phi[k];
    }
  }

}

template <class data_t>
ScalarField::Vars<data_t>::Vars()
{
  //Define the mapping from components of chombo grid to elements in Vars.
  //This allows to read/write data from the chombo grid into local
  //variables in Vars (which only exist for the current cell).

  //Scalars
  define_enum_mapping(c_chi, chi);
  define_enum_mapping(c_K, K);
  define_enum_mapping(c_Theta, Theta);
  define_enum_mapping(c_lapse, lapse);

  //Vectors
  define_enum_mapping(Interval(c_Gamma1,c_Gamma3), Gamma);
  define_enum_mapping(Interval(c_shift1,c_shift3), shift);
  define_enum_mapping(Interval(c_B1,c_B3), B);

  //Symmetric 2-tensors
  define_symmetric_enum_mapping(Interval(c_h11,c_h33), h);
  define_symmetric_enum_mapping(Interval(c_A11,c_A33), A);

  //Scalars - matter
  define_enum_mapping(c_phi, phi);   //Note that the matter fields are added here
  define_enum_mapping(c_Pi, Pi);

}

#endif /* SCALARFIELD_IMPL_HPP_ */
