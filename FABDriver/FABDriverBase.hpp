#ifndef FABDRIVERBASE_HPP_
#define FABDRIVERBASE_HPP_
//FABDriver and FABDriverBase implement looping of all points inside an FArrayBox
//with vectorisation and OpenMP.
//
//The reason why FABDriverBase is necessary is that if we want to inherit from one of the compute classes
//then the templating of the FABDriver with class compute_t becomes a problem when constructing the objects.
//FABDriverBase gives us all we need on the side of the computer but isn't templated.

//The following file has to be provided by the user.
//It must include c_NUM, the total number of simulation components.
#include "user_enum.hpp"


#include "tensor.hpp"

struct FABDriverBase
{
    const double *m_in_ptr[c_NUM];
    const int *m_in_lo;
    const int *m_in_hi;
    int m_stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];

    template <class data_t>
    void local_vars(int idx, data_t (&out)[c_NUM]) const;

    template <class data_t>
    void diff1(int idx, int stride, double dx, data_t (&out)[c_NUM]) const;

    template <class data_t>
    void diff2(int idx, int stride, double dx, data_t (&out)[c_NUM]) const;

    template <class data_t>
    void mixed_diff2(int idx, int stride1, int stride2, double dx, data_t (&out)[c_NUM]) const;

    template <class data_t>
    void advection(int idx, double dx, const tensor<1, data_t>& vec, data_t (&out)[c_NUM]) const;

    template <class data_t>
    void dissipation(int idx, double dx, data_t (&out)[c_NUM]) const;
};
#endif /* FABDRIVERBASE_HPP_ */
