#if !defined(FABDRIVER_HPP_)
#error "This file should only be included through FABDriver.hpp"
#endif

#ifndef FABDRIVER_IMPL_HPP_
#define FABDRIVER_IMPL_HPP_

#include "simd.hpp"

template <class compute_t>
template <typename... param_types>
FABDriver<compute_t>::FABDriver(param_types... params) : m_compute(*this, std::forward<param_types>(params)...) {}

template <class compute_t>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out, const Box& loop_box)
{
    //Note: this function can only be called if the compute class supports vectorisation.
    //i.e. it must be templated over the data type and allow double and simd<double> as data type.

    //Makes sure we are not requesting data outside the box of 'out'
    CH_assert(out.box().contains(loop_box));

    set_pointers(in,out);

    const int* loop_lo = loop_box.loVect();
    const int* loop_hi = loop_box.hiVect();

#pragma omp parallel for default(shared) collapse(CH_SPACEDIM-1)
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif
#if CH_SPACEDIM >= 3
    for (int z = loop_lo[2]; z <= loop_hi[2]; ++z)
#endif
        for (int y = loop_lo[1]; y <= loop_hi[1]; ++y)
        {
            int x_simd_max = loop_lo[0] + simd<double>::simd_len * (((loop_hi[0] - loop_lo[0] + 1) / simd<double>::simd_len) - 1);
            // SIMD LOOP
#pragma novector
            for (int x = loop_lo[0]; x <= x_simd_max; x += simd<double>::simd_len)
            {
                m_compute.template compute<simd<double> >(x,y,z);
            }
            // REMAINDER LOOP
#pragma novector
            for (int x = x_simd_max + simd<double>::simd_len; x <= loop_hi[0]; ++x)
            {
                m_compute.template compute<double>(x,y,z);
            }
        }
}

template <class compute_t>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out, const Box& loop_box, disable_simd)
{
    //Makes sure we are not requesting data outside the box of 'out'
    CH_assert(out.box().contains(loop_box));

    set_pointers(in,out);

    const int* loop_lo = loop_box.loVect();
    const int* loop_hi = loop_box.hiVect();

#pragma omp parallel for default(shared) collapse(CH_SPACEDIM-1)
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif
#if CH_SPACEDIM >= 3
    for (int z = loop_lo[2]; z <= loop_hi[2]; ++z)
#endif
        for (int y = loop_lo[1]; y <= loop_hi[1]; ++y)
        {
#pragma novector
            for (int x = loop_lo[0]; x <= loop_hi[0]; ++x)
            {
                m_compute.compute(x,y,z);
            }
        }
}

template <class compute_t>
void
FABDriver<compute_t>::set_pointers(const FArrayBox& in, FArrayBox& out)
{
    // dataPtr in Chombo does CH_assert bound check
    // which we don't want to do in a loop
    for (int i = 0; i < c_NUM; ++i)
    {
        m_in_ptr[i] = in.dataPtr(i);
        m_out_ptr[i] = out.dataPtr(i);
    }

    m_in_lo = in.loVect();
    m_in_hi = in.hiVect();
    m_in_stride[0] = 1;
    m_in_stride[1] = m_in_hi[0]-m_in_lo[0]+1;
#if CH_SPACEDIM >= 3
    m_in_stride[2] = (m_in_hi[1]-m_in_lo[1]+1)*m_in_stride[1];
#endif
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif

    m_out_lo = out.loVect();
    m_out_hi = out.hiVect();

    m_out_stride[0] = 1;
    m_out_stride[1] = m_out_hi[0]-m_out_lo[0]+1;
#if CH_SPACEDIM >= 3
    m_out_stride[2] = (m_out_hi[1]-m_out_lo[1]+1)*m_out_stride[1];
#endif
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif
}

template <class compute_t>
template <typename... simd_info>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out, simd_info... info)
{
    execute(in,out,out.box(), std::forward<simd_info>(info)...);
}

template <class compute_t>
template <typename... simd_info>
void
FABDriver<compute_t>::execute(const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fillGhosts, simd_info... info)
{
    DataIterator dit0  = in.dataIterator();
    int nbox = dit0.size();
    for(int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const FArrayBox& in_fab = in[di];
        FArrayBox& out_fab = out[di];

        Box out_box;
        if (fillGhosts) out_box = out_fab.box();
        else out_box = in.disjointBoxLayout()[di];

        execute(in_fab,out_fab,out_box, std::forward<simd_info>(info)...);
    }
}

#endif /* FABDRIVER_IMPL_HPP_ */
