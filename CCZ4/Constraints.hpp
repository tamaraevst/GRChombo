//This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriverBase.hpp"

#include "CCZ4Geometry.hpp"

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
        data_t Ham;
        tensor<1, data_t> Mom;
    };

protected:
    double m_dx;

public:
    Constraints(const FABDriverBase& driver, double dx);

    template <class data_t>
    void compute(int x, int y, int z);

protected:
    const FABDriverBase& m_driver;

    template <class data_t>
    void demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out);

    template <class data_t>
    void constraint_equations(vars_t<data_t> &vars,
             const vars_t<data_t> (&d1)[CH_SPACEDIM],
             const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]);
};

#include "Constraints.tpp"

#endif /* CONSTRAINTS_HPP_ */
