#ifndef FABDRIVER_HPP_
#define FABDRIVER_HPP_
//FABDriver and FABDriverBase implement looping of all points inside an FArrayBox
//with vectorisation and OpenMP.
//
//The reason why FABDriverBase is necessary is that if we want to inherit from one of the compute classes
//then the templating of the FABDriver with class compute_t becomes a problem when constructing the objects.
//FABDriverBase gives us all we need on the side of the computer but isn't templated.

#include "FABDriverBase.hpp"
#include "FArrayBox.H"
#include "tensor.hpp"

template <class compute_t>
class FABDriver : public FABDriverBase
{
public:
    compute_t m_compute;

    template <typename... param_types>
    FABDriver(param_types... params) : m_compute(compute_t(std::forward<param_types>(params)..., *this)) {};

    void execute(const FArrayBox& in, FArrayBox& out);

//    template <class data_t>
//    void local_vars(int idx, typename compute_t::template vars_t<data_t>& out);
//
//    template <class data_t>
//    void diff1(int idx, int stride, typename compute_t::template vars_t<data_t>& out);
//
//    template <class data_t>
//    void diff2(int idx, int stride, typename compute_t::template vars_t<data_t>& out);
//
//    template <class data_t>
//    void mixed_diff2(int idx, int stride1, int stride2, typename compute_t::template vars_t<data_t>& out);
//
//    template <class data_t>
//    void advection(int idx, const tensor<1, data_t>& vec, typename compute_t::template vars_t<data_t>& out);
//
//    template <class data_t>
//    void dissipation(int idx, typename compute_t::template vars_t<data_t>& out);
};

#include "FABDriver.tpp"

#endif /* FABDRIVER_HPP_ */
