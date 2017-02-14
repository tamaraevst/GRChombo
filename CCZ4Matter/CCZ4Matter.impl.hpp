// Last edited K Clough 31.01.17

#if !defined(CCZ4MATTER_HPP_)
#error "This file should only be included through CCZ4Matter.hpp"
#endif


#ifndef CCZ4MATTER_IMPL_HPP_
#define CCZ4MATTER_IMPL_HPP_

#define COVARIANTZ4

//TODO KClough: inline??
template <class matter_t>
CCZ4Matter<matter_t>::CCZ4Matter(const FABDriverBase& driver,
    params_t params,
    const typename matter_t::matter_params_t matter_params,
    double dx,
    double sigma,
    int formulation,
    double G_Newton)
    : CCZ4(driver, params, dx, sigma, formulation, 0.0), //No cosmological const
      m_matter_params (matter_params) , m_G_Newton (G_Newton) {}


template <class matter_t>
template <class data_t>
void CCZ4Matter<matter_t>::compute(int ix, int iy, int iz)
{
  //copy data from chombo gridpoint into local variables
  Vars<data_t> matter_vars;
  m_driver.local_vars(matter_vars);

  //work out first derivatives of variables on grid
  Vars< tensor<1, data_t> > d1;
  FOR1(idir) m_deriv.diff1(d1, idir);

  // Repeated derivatives
  // Work out second derivatives of variables on grid
  Vars< tensor<2,data_t> > d2;
  FOR1(idir) m_deriv.diff2(d2, idir);
  // Mixed derivatives
  // Note: no need to symmetrise explicitely, this is done in mixed_diff2
  m_deriv.mixed_diff2(d2, 1, 0);
  m_deriv.mixed_diff2(d2, 2, 0);
  m_deriv.mixed_diff2(d2, 2, 1);

  // Calculate advection terms
  Vars<data_t> advec;
  advec.assign(0.);
  FOR1(idir) m_deriv.add_advection(advec, matter_vars.shift[idir], idir);

  // Call CCZ4 RHS - work out RHS without matter, no dissipation
  Vars<data_t> matter_rhs;
  matter_rhs.assign(0.);
  matter_rhs = rhs_equation(matter_vars, d1, d2, advec);

  //add RHS matter terms from EM tensor for matter fields
  add_EMTensor_rhs(matter_rhs, matter_vars, d1, d2, advec);

  //TODO: K Clough
  //It is necessary to reassign the rhs to a new var here
  //but I do not know why, for now it doesn't cost much, so ok.
  Vars<data_t> total_rhs;
  total_rhs.assign(0.);

  total_rhs.chi = matter_rhs.chi;
  total_rhs.K = matter_rhs.K;
  total_rhs.Theta = matter_rhs.Theta;
  total_rhs.lapse = matter_rhs.lapse;

  FOR1(i)
  {
    total_rhs.Gamma[i] = matter_rhs.Gamma[i];
    total_rhs.shift[i] = matter_rhs.shift[i];
    total_rhs.B[i]     = matter_rhs.B[i];

    FOR1(j)
    {
      total_rhs.h[i][j]   =  matter_rhs.h[i][j];
      total_rhs.A[i][j]   =  matter_rhs.A[i][j];
    }
  }

  //add evolution of matter fields themselves
  matter_t my_matter(m_matter_params);
  my_matter.add_matter_rhs(total_rhs, matter_vars, d1, d2, advec);

  //Add dissipation to all terms
  FOR1(idir) m_deriv.add_dissipation(total_rhs, m_sigma, idir);

  //Write the rhs into the output FArrayBox
  m_driver.store_vars(total_rhs);

}

// Function to add in matter terms to CCZ4 rhs
template <class matter_t>
template <class data_t>
void CCZ4Matter<matter_t>::add_EMTensor_rhs(
    Vars<data_t>  &matter_rhs,
    const Vars<data_t>  &matter_vars,
    const Vars< tensor<1,data_t> >& d1,
    const Vars< tensor<2,data_t> >& d2,
    const Vars<data_t>  &advec) {

  using namespace TensorAlgebra;
 
  auto h_UU = compute_inverse(matter_vars.h);
  auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

  //Calculate elements of the decomposed stress energy tensor
  matter_t my_matter(m_matter_params);
  auto emtensor =  my_matter.compute_emtensor(matter_vars, d1, h_UU, chris.ULL, advec);

  // Update RHS for K and Theta depending on formulation
  if (m_formulation == USE_BSSN) {
    matter_rhs.K += 4.0*M_PI*m_G_Newton*matter_vars.lapse*(emtensor.S + emtensor.rho);
    matter_rhs.Theta += 0.0;
  } else {
    matter_rhs.K += 4.0*M_PI*m_G_Newton*matter_vars.lapse*(emtensor.S - 3*emtensor.rho);
    matter_rhs.Theta += - 8.0*M_PI*m_G_Newton*matter_vars.lapse*emtensor.rho;
  }

  // Update RHS for other variables
  tensor<2, data_t> Sij_TF = emtensor.Sij;
  make_trace_free(Sij_TF, matter_vars.h, h_UU);

  FOR2(i,j)
  {
    matter_rhs.A[i][j] +=
        - 8.0*M_PI*m_G_Newton*matter_vars.chi*matter_vars.lapse*Sij_TF[i][j];
  }

  FOR1(i)
  {
    data_t add_rhs_Gamma = 0.0;
    FOR1(j)
    {
      add_rhs_Gamma +=
          - 16.0*M_PI*m_G_Newton*matter_vars.lapse*h_UU[i][j]*emtensor.Si[j];
    }

    matter_rhs.Gamma[i] += add_rhs_Gamma;
    matter_rhs.B[i] += add_rhs_Gamma;

  }

}


#endif /* CCZ4MATTER_IMPL_HPP_ */
