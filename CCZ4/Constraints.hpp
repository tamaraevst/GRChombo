//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class Constraints
{
protected:
    template <class data_t>
    struct vars_t
    {
        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;

        vars_t(){}

        template <class arr_t>
        vars_t(const arr_t& in);
    };

    template <class data_t>
    struct constraints_t
    {
        data_t Ham;
        tensor<1, data_t> Mom;
    };

    const FABDriverBase& m_driver;
    const FourthOrderDerivatives m_deriv;

public:
    Constraints(const FABDriverBase& driver, double dx);

    template <class data_t>
    void compute(int x, int y, int z);

protected:
    template <class data_t>
    constraints_t<data_t>
    constraint_equations(
        vars_t<data_t> &vars,
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]
    );
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
