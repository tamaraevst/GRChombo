//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class Constraints
{
protected:
    template <class data_t>
    struct vars_t : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later
        using VarsBase<data_t>::define_symmetric_enum_mapping; //Saves us some writing later

        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;

        vars_t();
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
        const vars_t< tensor<1,data_t> >& d1,
        const vars_t< tensor<2,data_t> >& d2
    );
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
