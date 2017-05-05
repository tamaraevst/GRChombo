#if !defined(BOXLOOPS_HPP_)
#error "This file should only be included through BoxLoops.hpp"
#endif

#ifndef BOXLOOPS_IMPL_HPP_
#define BOXLOOPS_IMPL_HPP_

#include "simd.hpp"
#include "Cell.hpp"
#include "BoxPointers.hpp"
#include "DebuggingTools.hpp"

template <typename... compute_ts>
ALWAYS_INLINE
void
BoxLoops::innermost_loop(ComputeClassPack<compute_ts...> compute_class_pack, const BoxPointers& box_pointers,
                         const int iy, const int iz, const int loop_lo_x, const int loop_hi_x)
{
#ifdef EQUATION_DEBUG_MODE /*TODO: MK: think of a nicer way. This macro jungle is a bit too much.*/
    innermost_loop(compute_class, iy, iz, loop_lo_x, loop_hi_x, disable_simd()); //In equation debug mode never use simd
#else
    int x_simd_max = loop_lo_x + simd<double>::simd_len * (((loop_hi_x - loop_lo_x + 1) / simd<double>::simd_len) - 1);
    // SIMD LOOP
#pragma novector
    for (int ix = loop_lo_x; ix <= x_simd_max; ix += simd<double>::simd_len)
    {
        //If you were sent here by a compile error of no matching function call make sure that
        //the compute class you are using allows for vectorisation (is templated over the data type)
        //To switch vectorisation off in a vectorised compute class pass disable_simd as last parameter to the
        //loop function. For a compute class without simd support pass no_simd_support().
        compute_class_pack.call_compute(Cell(IntVect(ix,iy,iz), box_pointers) );
    }
    // REMAINDER LOOP
#pragma novector
    for (int ix = x_simd_max + simd<double>::simd_len; ix <= loop_hi_x; ++ix)
    {
        compute_class_pack.call_compute(Cell(IntVect(ix,iy,iz), box_pointers), disable_simd());
    }
#endif
}

template <typename... compute_ts, typename... simd_info>
ALWAYS_INLINE
void
BoxLoops::innermost_loop(ComputeClassPack<compute_ts...> compute_class_pack, const BoxPointers& box_pointers,
                         const int iy, const int iz, const int loop_lo_x, const int loop_hi_x, simd_info... info)
{
#pragma novector
    for (int ix = loop_lo_x; ix <= loop_hi_x; ++ix)
    {
#ifdef EQUATION_DEBUG_MODE
        EquationDebugging::set_global_cell_coordinates(IntVect(ix,iy,iz));
#endif
        compute_class_pack.call_compute( Cell(IntVect(ix,iy,iz), box_pointers), std::forward<simd_info>(info)... );
    }
}

template <typename... compute_ts, typename... simd_info>
void
BoxLoops::loop(ComputeClassPack<compute_ts...> compute_class_pack, const FArrayBox& in, FArrayBox& out, const Box& loop_box, simd_info... info)
{
    //Makes sure we are not requesting data outside the box of 'out'
    CH_assert(out.box().contains(loop_box));

    BoxPointers box_pointers(in,out);

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
            innermost_loop(compute_class_pack, box_pointers, iy, iz, loop_lo[0], loop_hi[0], std::forward<simd_info>(info)...);
        }
}

template <typename compute_t, typename... simd_info>
void
BoxLoops::loop(compute_t compute_class, const FArrayBox& in, FArrayBox& out, const Box & loop_box, simd_info... info)
{
    loop(make_compute_pack(compute_class), in, out, loop_box, std::forward<simd_info>(info)...);
}

template <typename... compute_ts, typename... simd_info>
void
BoxLoops::loop(ComputeClassPack<compute_ts...> compute_class_pack, const FArrayBox& in, FArrayBox& out, simd_info... info)
{
    loop(compute_class_pack, in,out,out.box(), std::forward<simd_info>(info)...);
}

template <typename compute_t, typename... simd_info>
void
BoxLoops::loop(compute_t compute_class, const FArrayBox& in, FArrayBox& out, simd_info... info)
{
    loop(make_compute_pack(compute_class), in, out, std::forward<simd_info>(info)...);
}

template <typename... compute_ts, typename... simd_info>
void
BoxLoops::loop(ComputeClassPack<compute_ts...> compute_class_pack, const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fill_ghosts, simd_info... info)
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

        loop(compute_class_pack, in_fab,out_fab,out_box, std::forward<simd_info>(info)...);
    }
}

template <typename compute_t, typename... simd_info>
void
BoxLoops::loop(compute_t compute_class, const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fill_ghosts, simd_info... info)
{
    loop(make_compute_pack(compute_class), in, out, fill_ghosts, std::forward<simd_info>(info)...);
}

#endif /* BOXLOOPS_IMPL_HPP_ */
