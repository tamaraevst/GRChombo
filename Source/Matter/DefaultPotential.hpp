// Last edited K Clough 13.05.17

#ifndef DEFAULTPOTENTIAL_HPP_
#define DEFAULTPOTENTIAL_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include <array>

class DefaultPotential {
public:
    //! The constructor
    DefaultPotential() {}

    //! Set the potential function for the scalar field here to zero
    template <class data_t>
    void
    compute_potential(data_t &V_of_phi, data_t &dVdphi, const data_t phi_here) const
    {
        //The potential value at phi
        V_of_phi = 0.0;

        //The potential gradient at phi
        dVdphi = 0.0;
    }
};

#endif /* DEFAULTPOTENTIAL_HPP_ */
