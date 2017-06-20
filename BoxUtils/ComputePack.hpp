#ifndef COMPUTEPACK_HPP_
#define COMPUTEPACK_HPP_

#include "DebuggingTools.hpp"
#include "Cell.hpp"

template <typename... compute_ts>
class ComputePack
{
    std::tuple<compute_ts...> m_compute_tuple;

    //Begin: Helper functions for calling 'compute' for several compute classes -->
    template<std::size_t ID = 0, typename... helper_compute_ts, class data_t>
    ALWAYS_INLINE
    typename std::enable_if<ID == sizeof...(helper_compute_ts), void>::type
    call_compute_helper(std::tuple<helper_compute_ts...>&, Cell<data_t> current_cell) { } //If we have reached the end of the tuple do nothing

    template<std::size_t ID = 0, typename... helper_compute_ts, class data_t>
    ALWAYS_INLINE
    typename std::enable_if<ID < sizeof...(helper_compute_ts), void>::type
    call_compute_helper(std::tuple<helper_compute_ts...>& compute_pack, Cell<data_t> current_cell)
    {
        std::get<ID>(compute_pack).compute(current_cell); //Call compute for the current component
        call_compute_helper<ID + 1, helper_compute_ts...>(compute_pack, current_cell); //call again for next component
    }
    //End: Helper functions for calling 'compute' for several compute classes

public:
    ComputePack(const std::tuple<compute_ts...>& compute_tuple) : m_compute_tuple (compute_tuple) {}

    template <class data_t, typename... simd_info>
    void call_compute(const Cell<data_t>& current_cell, simd_info... info)
    {
        call_compute_helper(m_compute_tuple, current_cell, std::forward<simd_info>(info)...);
    }
};

template <typename... compute_ts>
ComputePack<compute_ts...> make_compute_pack(compute_ts... compute_classes)
{
    return ComputePack<compute_ts...>(std::make_tuple(std::forward<compute_ts>(compute_classes)...));
}

#endif /* COMPUTEPACK_HPP_ */
