//This class enforces A to be trace-free
#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "UserVariables.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"
#include "Interval.H"

#include <array>

class EnforceTfA
{
protected:
    const FABDriverBase& m_driver;

public:
    EnforceTfA(const FABDriverBase& driver) :
        m_driver (driver)
    {}

    template <class data_t>
    struct Vars : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_symmetric_enum_mapping; //Saves us some writing later
        Vars();

        tensor<2, data_t> h;
        tensor<2, data_t> A;
    };

    template <class data_t>
    void compute(int ix, int iy, int iz)
    {
        Vars<data_t> vars;
        m_driver.local_vars(vars);

        auto h_UU = TensorAlgebra::compute_inverse(vars.h);
        TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

        m_driver.store_vars(vars, Interval(c_A11, c_A33));
    }
};


template <class data_t>
EnforceTfA::Vars<data_t>::Vars()
{
    //Symmetric 2-tensors
    define_symmetric_enum_mapping(Interval(c_h11,c_h33), h);
    define_symmetric_enum_mapping(Interval(c_A11,c_A33), A);
}

#endif /* FIXTFA_HPP_ */
