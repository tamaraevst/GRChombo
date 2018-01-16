#ifndef COMPUTEPACK_HPP_
#define COMPUTEPACK_HPP_

#include "DebuggingTools.hpp"
#include "Cell.hpp"

template <typename... compute_ts>
class ComputePack
{
    std::tuple<compute_ts...> m_compute_tuple;

    //Begin: Helper functions for calling 'compute' for several compute classes -->
    template<std::size_t ID = 0, class data_t>
    ALWAYS_INLINE
    typename std::enable_if<ID == sizeof...(compute_ts), void>::type
    call_compute_helper(Cell<data_t> current_cell) const { } //If we have reached the end of the tuple do nothing

    template<std::size_t ID = 0, class data_t>
    ALWAYS_INLINE
    typename std::enable_if<ID < sizeof...(compute_ts), void>::type
    call_compute_helper(Cell<data_t> current_cell) const
    {
        std::get<ID>(m_compute_tuple).compute(current_cell); //Call compute for the current component
        call_compute_helper<ID + 1>(current_cell); //call again for next component
    }
    //End: Helper functions for calling 'compute' for several compute classes

public:
    ComputePack(const std::tuple<compute_ts...>& compute_tuple) : m_compute_tuple (compute_tuple) {}

    template <class data_t>
    void call_compute(const Cell<data_t>& current_cell) const
    {
        call_compute_helper(current_cell);
    }
};

///Helper struct for checking whether a given template argument is a compute pack
template <typename... Ts>
struct is_compute_pack : public std::false_type {};

template <typename... Ts>
struct is_compute_pack<ComputePack<Ts...>> : public std::true_type {};
//End: Helper struct for checking whether a given template argument is a compute pack

template <typename... compute_ts>
ComputePack<compute_ts...> make_compute_pack(compute_ts... compute_classes)
{
    return ComputePack<compute_ts...>(std::make_tuple(std::forward<compute_ts>(compute_classes)...));
}

#endif /* COMPUTEPACK_HPP_ */
