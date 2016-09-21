//This class enforces A to be trace-free
#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "UserVariables.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"

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
    struct vars_t : VarsBase<data_t>
    {
        using VarsBase<data_t>::m_assignment_ptrs; //Saves us some writing later

        tensor<2, data_t> h;
        tensor<2, data_t> A;

        vars_t()
        {
            m_assignment_ptrs[c_h11].push_back(&h[0][0]);
            m_assignment_ptrs[c_h12].push_back(&h[0][1]);
            m_assignment_ptrs[c_h12].push_back(&h[1][0]);
            m_assignment_ptrs[c_h13].push_back(&h[0][2]);
            m_assignment_ptrs[c_h13].push_back(&h[2][0]);
            m_assignment_ptrs[c_h22].push_back(&h[1][1]);
            m_assignment_ptrs[c_h23].push_back(&h[1][2]);
            m_assignment_ptrs[c_h23].push_back(&h[2][1]);
            m_assignment_ptrs[c_h33].push_back(&h[2][2]);

            m_assignment_ptrs[c_A11].push_back(&A[0][0]);
            m_assignment_ptrs[c_A12].push_back(&A[0][1]);
            m_assignment_ptrs[c_A12].push_back(&A[1][0]);
            m_assignment_ptrs[c_A13].push_back(&A[0][2]);
            m_assignment_ptrs[c_A13].push_back(&A[2][0]);
            m_assignment_ptrs[c_A22].push_back(&A[1][1]);
            m_assignment_ptrs[c_A23].push_back(&A[1][2]);
            m_assignment_ptrs[c_A23].push_back(&A[2][1]);
            m_assignment_ptrs[c_A33].push_back(&A[2][2]);
        };
    };

    template <class data_t>
    void compute(int x, int y, int z)
    {
        idx_t<data_t> idx = m_driver.in_idx(x, y, z);

        vars_t<data_t> vars;
        m_driver.local_vars(vars,idx);

        auto h_UU = TensorAlgebra::compute_inverse(vars.h);
        TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

        idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A11])[out_idx] = vars.A[0][0];
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A12])[out_idx] = vars.A[0][1];
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A13])[out_idx] = vars.A[0][2];
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A22])[out_idx] = vars.A[1][1];
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A23])[out_idx] = vars.A[1][2];
        SIMDIFY<data_t>(m_driver.m_out_ptr[c_A33])[out_idx] = vars.A[2][2];
    }
};

#endif /* FIXTFA_HPP_ */
