#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriver.hpp"

class FixTfA
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
    FixTfA(double dx, FABDriver& driver);

protected:
    template <class data_t>
    void compute(int x, int y, int z);

    template <class data_t>
    void fix_tf_A(vars_t<data_t> &vars)
    {
       auto h_UU = compute_inverse_metric(vars);
       make_trace_free(vars.A, vars.h, h_UU);
    }
};

class CCZ4Factory
{


}

#endif /* FIXTFA_HPP_ */
