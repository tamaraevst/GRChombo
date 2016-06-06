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
#include "LevelData.H"
#include "tensor.hpp"

template <class compute_t>
class FABDriver : public FABDriverBase
{
public:
    compute_t m_compute;

    template <typename... param_types>
    FABDriver(param_types... params) :
        m_compute(std::forward<param_types>(params)..., *this)
    {};

    void execute(const FArrayBox& in, FArrayBox& out);

    void execute(const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out);
};

#include "FABDriver.tpp"

#endif /* FABDRIVER_HPP_ */
