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
    FABDriver(param_types... params);

    //Takes input 'in', writes output into the subox 'loop_box' of 'out'
    void execute(const FArrayBox& in, FArrayBox& out, const Box & loop_box);
    void execute(const FArrayBox& in, FArrayBox& out, const Box & loop_box, disable_simd);

    template <typename... simd_info>
    void execute(const FArrayBox& in, FArrayBox& out, simd_info... info); //Uses out.box() as loop_box

    //MK: Could give the ghost treatment a default argument but I think it's better to force the user to make a concious decision
    //Wrong fillGhosts can give errors that are very hard to debug
    template <typename... simd_info>
    void execute(const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fillGhosts, simd_info... info);

    void set_pointers(const FArrayBox& in, FArrayBox& out);
};

#include "FABDriver.impl.hpp"

#endif /* FABDRIVER_HPP_ */
