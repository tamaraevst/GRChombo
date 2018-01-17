// Last edited K Clough 09.05.17

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include <array>

class Potential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           data_t phi_here) const
    {
        // The potential value at phi
        V_of_phi = pow(m_params.scalar_mass * phi_here, 2.0); // e.g. m^2 phi^2

        // The potential gradient at phi
        dVdphi =
            2.0 * pow(m_params.scalar_mass, 2.0) * phi_here; //  e.g. 2 m^2 phi
    }
};

#endif /* POTENTIAL_HPP_ */
