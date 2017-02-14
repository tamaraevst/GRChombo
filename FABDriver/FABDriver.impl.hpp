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
ALWAYS_INLINE
void
FABDriver<compute_t>::innermost_loop(const int iy, const int iz, const int loop_lo_x, const int loop_hi_x)
{
    int x_simd_max = loop_lo_x + simd<double>::simd_len * (((loop_hi_x - loop_lo_x + 1) / simd<double>::simd_len) - 1);
    // SIMD LOOP
#pragma novector
    for (int ix = loop_lo_x; ix <= x_simd_max; ix += simd<double>::simd_len)
    {
        //If you were sent here by a compile error of no matching function call make sure that
        //the compute class you are using allows for vectorisation (is templated over the data type)
        //To switch vectorisation off in a vectorised compute class pass disable_simd as last parameter to the
        //execute function. For a compute class without simd support pass no_simd_support().
        set_idx(ix,iy,iz);
        m_compute.template compute<simd<double> >(ix,iy,iz);
    }
    // REMAINDER LOOP
#pragma novector
    for (int ix = x_simd_max + simd<double>::simd_len; ix <= loop_hi_x; ++ix)
    {
        set_idx(ix,iy,iz);
        m_compute.template compute<double>(ix,iy,iz);
    }
}

template <class compute_t>
ALWAYS_INLINE
void
FABDriver<compute_t>::innermost_loop(const int iy, const int iz, const int loop_lo_x, const int loop_hi_x, disable_simd)
{
#pragma novector
    for (int ix = loop_lo_x; ix <= loop_hi_x; ++ix)
    {
        set_idx(ix,iy,iz);
        m_compute.template compute<double>(ix,iy,iz);
    }
}

template <class compute_t>
ALWAYS_INLINE
void
FABDriver<compute_t>::innermost_loop(const int iy, const int iz, const int loop_lo_x, const int loop_hi_x, no_simd_support)
{
#pragma novector
    for (int ix = loop_lo_x; ix <= loop_hi_x; ++ix)
    {
        set_idx(ix,iy,iz);
        m_compute.compute(ix,iy,iz);
    }
}

template <class compute_t>
template <typename... simd_info>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out, const Box& loop_box, simd_info... info)
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
    for (int iz = loop_lo[2]; iz <= loop_hi[2]; ++iz)
#endif
        for (int iy = loop_lo[1]; iy <= loop_hi[1]; ++iy)
        {
            innermost_loop(iy, iz, loop_lo[0], loop_hi[0], std::forward<simd_info>(info)...);
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
FABDriver<compute_t>::execute(const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fill_ghosts, simd_info... info)
{
    DataIterator dit0  = in.dataIterator();
    int nbox = dit0.size();
    for(int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const FArrayBox& in_fab = in[di];
        FArrayBox& out_fab = out[di];

        Box out_box;
        if (fill_ghosts) out_box = out_fab.box();
        else out_box = in.disjointBoxLayout()[di];

        execute(in_fab,out_fab,out_box, std::forward<simd_info>(info)...);
    }
}

#endif /* FABDRIVER_IMPL_HPP_ */
