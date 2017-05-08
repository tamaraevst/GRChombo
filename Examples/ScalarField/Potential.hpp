// Last edited K Clough 31.01.17

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"
#include "tensor.hpp"

#include <array>

//! Set the potential function for the scalar field here
template <class data_t>
void compute_potential(data_t &V_of_phi, data_t &dVdphi, const data_t phi_here, const double scalar_mass) {

  //The potential value at phi
  V_of_phi = pow(scalar_mass*phi_here, 2.0); // e.g. m^2 phi^2

  //The potential gradient at phi
  dVdphi = 2.0*pow(scalar_mass,2.0)*phi_here;  //  e.g. 2 m^2 phi

}

#endif /* POTENTIAL_HPP_ */
