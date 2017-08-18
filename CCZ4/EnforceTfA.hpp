//This class enforces A to be trace-free
#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "UserVariables.hpp"
#include "tensor.hpp"
#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"
#include "Interval.H"
#include "Cell.hpp"
#include "VarsTools.hpp"

#include <array>

class EnforceTfA
{
public:
    template <class data_t>
    struct Vars
    {
        tensor<2, data_t> h;
        tensor<2, data_t> A;

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function);
    };

    template <class data_t>
    void compute(Cell<data_t> current_cell)
    {
        auto vars = current_cell.template load_vars<Vars>();

        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
        TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

        current_cell.store_vars(vars);
    }
};

template <class data_t>
template <typename mapping_function_t>
void EnforceTfA::Vars<data_t>::enum_mapping(mapping_function_t mapping_function)
{
    VarsTools::define_symmetric_enum_mapping(mapping_function, GRInterval<c_h11,c_h33>(), h);
    VarsTools::define_symmetric_enum_mapping(mapping_function, GRInterval<c_A11,c_A33>(), A);
}

#endif /* FIXTFA_HPP_ */
