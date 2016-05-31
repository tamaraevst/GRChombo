#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriver.hpp"

class Constraints
{
protected:
    template <class data_t>
    struct vars_t   //TODO: All this needs to go elsewhere!
    {
        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;
    };

protected:
    double m_dx;

    const double *m_in_ptr[c_NUM];
    const int *m_in_lo;
    const int *m_in_hi;
    int m_stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];
    
public:
    Constraints(params_t params, double dx, FABDriver& driver);

protected:
    template <class data_t>
    void compute(int x, int y, int z);

    template <class data_t>
    void constraint_equations(const vars_t<data_t> &vars,
             const vars_t<data_t> (&d1)[CH_SPACEDIM],
             const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
             const vars_t<data_t> &advec,
             vars_t<data_t> &rhs);
};

class CCZ4Factory
{


}

#endif /* CONSTRAINTS_HPP_ */
