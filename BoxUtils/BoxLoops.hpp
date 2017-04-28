#ifndef BOXLOOPS_HPP_
#define BOXLOOPS_HPP_

#include "FArrayBox.H"
#include "LevelData.H"
#include "tensor.hpp"
#include "BoxPointers.hpp"

enum
{
    SKIP_GHOST_CELLS,
    FILL_GHOST_CELLS
};

namespace BoxLoops
{
    //Begin: Helper functions for calling 'compute' for several compute classes -->
    template<std::size_t ID = 0, typename... compute_ts>
    ALWAYS_INLINE
    typename std::enable_if<ID == sizeof...(compute_ts), void>::type
    call_compute(std::tuple<compute_ts...>&, Cell current_cell) { } //If we have reached the end of the tuple do nothing

    template<std::size_t ID = 0, typename... compute_ts>
    ALWAYS_INLINE
    typename std::enable_if<ID < sizeof...(compute_ts), void>::type
    call_compute(std::tuple<compute_ts...>& compute_pack, Cell current_cell)
    {
        std::get<ID>(compute_pack).template compute<simd<double>>(current_cell); //Call compute for the current component
        call_compute<ID + 1, compute_ts...>(compute_pack, current_cell); //call again for next component
    }
    //End: Helper functions for calling 'compute' for several compute classes

    ///Perform the innermost (i.e. the x) loop, standard version with simd
    template <typename... compute_ts>
    void innermost_loop(std::tuple<compute_ts...> compute_class_pack, const BoxPointers& box_pointers,
                        const int iy, const int iz, const int loop_lo_x, const int loop_hi_x);

    ///Perform the innermost (i.e. the x) loop, switches simd off for compute classes which have simd support
    template <typename... compute_ts>
    void innermost_loop(std::tuple<compute_ts...> compute_class_pack, const BoxPointers& box_pointers,
                        const int iy, const int iz, const int loop_lo_x, const int loop_hi_x, disable_simd);

    ///Perform the innermost (i.e. the x) loop, for compute classes without simd support
    ///(i.e. where the compute function is not templated at all.
    template <typename... compute_ts>
    void innermost_loop(std::tuple<compute_ts...> compute_class_pack, const BoxPointers& box_pointers,
                        const int iy, const int iz, const int loop_lo_x, const int loop_hi_x, no_simd_support);

    ///Performs loop insde the box loop_box and calls compute(...) for all compute classes in the compute_class_pack
    ///with input FArrayBox 'in' and output FArrayBox 'out'.
    template <typename... compute_ts, typename... simd_info>
    void loop(std::tuple<compute_ts...> compute_class_pack, const FArrayBox& in, FArrayBox& out, const Box & loop_box, simd_info... info);

    ///Same as above but for only one compute class (rather than a pack of them)
    template <typename compute_t, typename... simd_info>
    void loop(compute_t compute_class, const FArrayBox& in, FArrayBox& out, const Box & loop_box, simd_info... info);

    ///Performs loop insde the whole box of 'out' and calls compute(...) for all compute classes in the compute_class_pack
    ///with input FArrayBox 'in' and output FArrayBox 'out'.
    template <typename... compute_ts, typename... simd_info>
    void loop(std::tuple<compute_ts...> compute_class_pack, const FArrayBox& in, FArrayBox& out, simd_info... info); //Uses out.box() as loop_box

    ///Same as above but for only one compute class (rather than a pack of them)
    template <typename compute_t, typename... simd_info>
    void loop(compute_t compute_class, const FArrayBox& in, FArrayBox& out, simd_info... info); //Uses out.box() as loop_box

    ///Performs loop over all boxes and inside all boxes of the LevelData 'out' and calls compute(...) for all compute
    //classes in the compute_class_pack with input data taken from 'in' and output written to 'out'
    //MK: Could give the ghost treatment a default argument but I think it's better to force the user to make a concious decision
    //Wrong fill_ghosts can give errors that are very hard to debug
    template <typename... compute_ts, typename... simd_info>
    void loop(std::tuple<compute_ts...> compute_class_pack, const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fill_ghosts, simd_info... info);

    ///Same as above but for only one compute class (rather than a pack of them)
    template <typename compute_t, typename... simd_info>
    void loop(compute_t compute_class, const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fill_ghosts, simd_info... info);
}

#include "BoxLoops.impl.hpp"

#endif /* BOXLOOPS_HPP_ */
