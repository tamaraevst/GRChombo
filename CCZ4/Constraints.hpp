#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriver.hpp"

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
    };

protected:
    double m_dx;

public:
    Constraints(double dx, const FABDriverBase& driver);

protected:
    template <class data_t>
    void compute(int x, int y, int z);

    const FABDriverBase& m_driver;

    template <class data_t>
    void constraint_equations(vars_t<data_t> &vars,
             const vars_t<data_t> (&d1)[CH_SPACEDIM],
             const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM]);
};
#endif /* CONSTRAINTS_HPP_ */
