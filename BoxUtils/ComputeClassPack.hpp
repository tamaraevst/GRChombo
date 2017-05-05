#ifndef COMPUTECLASSPACK_HPP_
#define COMPUTECLASSPACK_HPP_

template <typename... compute_ts>
class ComputeClassPack
{
    std::tuple<compute_ts...> m_compute_tuple;

    //Begin: Helper functions for calling 'compute' for several compute classes -->
    template<std::size_t ID = 0, typename... helper_compute_ts, typename... forwarded_params_t>
    ALWAYS_INLINE
    typename std::enable_if<ID == sizeof...(helper_compute_ts), void>::type
    call_compute_helper(std::tuple<helper_compute_ts...>&, forwarded_params_t... forwarded_params) { } //If we have reached the end of the tuple do nothing

    template<std::size_t ID = 0, typename... helper_compute_ts, typename... forwarded_params_t>
    ALWAYS_INLINE
    typename std::enable_if<ID < sizeof...(helper_compute_ts), void>::type
    call_compute_helper(std::tuple<helper_compute_ts...>& compute_pack, forwarded_params_t... forwarded_params)
    {
        call_single_compute(std::get<ID>(compute_pack), std::forward<forwarded_params_t>(forwarded_params)...); //Call compute for the current component
        call_compute_helper<ID + 1, helper_compute_ts...>(compute_pack, std::forward<forwarded_params_t>(forwarded_params)...); //call again for next component
    }

    template <class compute_t>
    void call_single_compute(compute_t& compute_class, const Cell& current_cell)
    {
        compute_class.template compute<simd<double>>(current_cell);
    }

    template <class compute_t>
    void call_single_compute(compute_t& compute_class, const Cell& current_cell, disable_simd)
    {
        compute_class.template compute<double>(current_cell);
    }

    template <class compute_t>
    void call_single_compute(compute_t& compute_class, const Cell& current_cell, no_simd_support)
    {
        compute_class.compute(current_cell);
    }
    //End: Helper functions for calling 'compute' for several compute classes

public:
    ComputeClassPack(const std::tuple<compute_ts...>& compute_tuple) : m_compute_tuple (compute_tuple) {}

    template <typename... simd_info>
    void call_compute(const Cell& current_cell, simd_info... info)
    {
        call_compute_helper(m_compute_tuple, current_cell, std::forward<simd_info>(info)...);
    }
};

template <typename... compute_ts>
ComputeClassPack<compute_ts...> make_compute_pack(compute_ts... compute_classes)
{
    return ComputeClassPack<compute_ts...>(std::make_tuple(std::forward<compute_ts>(compute_classes)...));
}

#endif /* COMPUTECLASSPACK_HPP_ */
