#if !defined(SFMATTER_HPP_)
#error "This file should only be included through SFMatter.hpp"
#endif

#ifndef SFMATTER_IMPL_HPP_
#define SFMATTER_IMPL_HPP_

#define COVARIANTZ4

//inline

template <class data_t>
auto SFMatter::calc_emtensor(
    const vars_t<data_t> &vars,
    const vars_t< tensor<1,data_t> >& d1,
    const tensor<2, data_t> &h_UU,
    const tensor<3, data_t> &chris,
    const vars_t<data_t> &advec) -> emtensor_t<data_t> {

  emtensor_t<data_t> out;

  // Calculate the stress energy tensor elements

  // Find the potential and its gradient in terms of phi
  potential_t<data_t> potential = calc_potential(vars.phi);// vars.phi*vars.phi;

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
  FOR2(i,j)
  {
    out.Sij[i][j] = -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
  }

  out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij , h_UU);

  FOR1(i)
  {
    out.Si[i] = - T_i[i]/vars.lapse;

    FOR1(j)
    {
      out.Si[i] += vars.shift[j]/vars.lapse * out.Sij[i][j];
    }
  }

  out.rho = dphidt2/vars.lapse/vars.lapse + 0.5*Vt; // = T_ab * n^a n^b
  FOR2(i,j)
  {
    out.rho += (-0.5*Vt*vars.h[i][j]/vars.chi + out.Sij[i][j])
              * vars.shift[i]*vars.shift[j]/vars.lapse/vars.lapse;
  }

  return out;
}

/// Set the potential function for the scalar field here
template <class data_t>
auto SFMatter::calc_potential(const data_t phi) -> potential_t<data_t> {

  potential_t<data_t> out;

  out.V_of_phi = 0*phi*phi; // e.g. m^2 phi^2 NB:WOULD LIKE COSINES HERE
  out.dVdphi = 0*2*phi;  //  e.g. 2 m^2 phi

  return out;
}

template <class data_t>
auto SFMatter::calc_total_rhs(
    const CCZ4::vars_t<data_t> &CCZ4_rhs,
    const vars_t<data_t> &matter_rhs,
    const vars_t<data_t> &vars,
    const vars_t< tensor<1,data_t> >& d1,
    const vars_t< tensor<2,data_t> >& d2,
    const vars_t<data_t> &advec) -> vars_t<data_t> {

  //Set all RHS to zero, then add the CCZ4 non zero terms
  vars_t<data_t> total_rhs;
  total_rhs.assign(0.);

  // Need to add CCZ4 terms separately as data structure excludes the matter terms
  total_rhs.chi   += CCZ4_rhs.chi    + matter_rhs.chi;
  total_rhs.K     += CCZ4_rhs.K      + matter_rhs.K;
  total_rhs.Theta += CCZ4_rhs.Theta  + matter_rhs.Theta;
  total_rhs.lapse += CCZ4_rhs.lapse  + matter_rhs.lapse;

  FOR1(i)
  {
    total_rhs.Gamma[i] += CCZ4_rhs.Gamma[i] + matter_rhs.Gamma[i];
    total_rhs.shift[i] += CCZ4_rhs.shift[i] + matter_rhs.shift[i];
    total_rhs.B[i]     += CCZ4_rhs.B[i]     + matter_rhs.B[i];

    FOR1(j)
    {
      total_rhs.h[i][j]   += CCZ4_rhs.h[i][j] + matter_rhs.h[i][j];
      total_rhs.A[i][j]   += CCZ4_rhs.A[i][j] + matter_rhs.A[i][j];
    }
  }

  using namespace TensorAlgebra;

  auto h_UU = compute_inverse(vars.h);
  auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

  emtensor_t<data_t> emtensor = calc_emtensor(vars, d1, h_UU, chris.ULL, advec);
  potential_t<data_t> potential = calc_potential(vars.phi);

  //evolution equations for scalar field and (minus) its conjugate momentum
  total_rhs.phi = vars.lapse * vars.Pi + advec.phi;

  total_rhs.Pi = vars.lapse*(vars.K * vars.Pi - potential.dVdphi) + advec.Pi;

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

  return total_rhs;
}

template <class data_t>
SFMatter::vars_t<data_t>::vars_t()
{
  //Define the mapping from components of chombo grid to elements in vars_t.
  //This allows to read/write data from the chombo grid into local
  //variables in vars_t (which only exist for the current cell).

  //Scalars
  define_enum_mapping(c_chi, chi);
  define_enum_mapping(c_K, K);
  define_enum_mapping(c_Theta, Theta);
  define_enum_mapping(c_lapse, lapse);
  define_enum_mapping(c_phi, phi);   //Note that the matter fields are added here
  define_enum_mapping(c_Pi, Pi);

  //Vectors
  define_enum_mapping(Interval(c_Gamma1,c_Gamma3), Gamma);
  define_enum_mapping(Interval(c_shift1,c_shift3), shift);
  define_enum_mapping(Interval(c_B1,c_B3), B);

  //Symmetric 2-tensors
  define_symmetric_enum_mapping(Interval(c_h11,c_h33), h);
  define_symmetric_enum_mapping(Interval(c_A11,c_A33), A);
}

#endif /* SFMATTER_IMPL_HPP_ */
