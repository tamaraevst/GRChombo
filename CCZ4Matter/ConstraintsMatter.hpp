//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTSMATTER_HPP_
#define CONSTRAINTSMATTER_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"

#include "CCZ4Geometry.hpp"
#include "CCZ4EMTensorSF.hpp"

#include <array>

class ConstraintsMatter
{
protected:
    template <class data_t>
    struct vars_t : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later

        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        tensor<1, data_t> shift;
				data_t lapse;
        data_t phi;
        data_t PiM;

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
    ConstraintsMatter(const FABDriverBase& driver, double dx);

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

#include "ConstraintsMatter.impl.hpp"

#endif /* CONSTRAINTSMATTER_HPP_ */
