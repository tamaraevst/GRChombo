//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Cell.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class Constraints
{
protected:
    template <class data_t>
    struct Vars
    {
        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;

        /// Defines the mapping between members of Vars and Chombo grid variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t>
    struct constraints_t
    {
        data_t Ham;
        tensor<1, data_t> Mom;
    };

    const FourthOrderDerivatives m_deriv;
    double m_cosmological_constant;

public:
    Constraints(double dx, double cosmological_constant = 0);

    template <class data_t>
    void compute(Cell<data_t> current_cell);

protected:
    template <class data_t, template<typename> class vars_t, template<typename> class diff2_vars_t>
    constraints_t<data_t>
    constraint_equations(
        const vars_t<data_t> &vars,
        const vars_t< tensor<1,data_t> >& d1,
        const diff2_vars_t< tensor<2,data_t> >& d2
    );
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
