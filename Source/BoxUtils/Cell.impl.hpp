#if !defined(CELL_HPP_)
#error "This file should only be included through Cell.hpp"
#endif

#ifndef CELL_IMPL_HPP_
#define CELL_IMPL_HPP_

#include "GRInterval.hpp"
#include "simd.hpp"

template <class data_t>
ALWAYS_INLINE data_t Cell<data_t>::load_vars(const int ivar) const
{
    return SIMDIFY<data_t>(m_box_pointers.m_in_ptr[ivar])[m_in_index];
}

template <class data_t>
ALWAYS_INLINE void Cell<data_t>::load_vars(data_t &out, const int ivar) const
{
    out = load_vars(ivar);
}

template <class data_t>
void Cell<data_t>::load_vars(data_t (&out)[NUM_VARS]) const
{
    for (int ivar = 0; ivar < NUM_VARS; ++ivar)
    {
        out[ivar] = load_vars(ivar);
    }
}

template <class data_t>
template <template <typename> class vars_t>
void Cell<data_t>::load_vars(vars_t<data_t> &vars) const
{
    vars.enum_mapping([&](const int &ivar, data_t &var) {
        var = SIMDIFY<data_t>(m_box_pointers.m_in_ptr[ivar])[m_in_index];
    });
}

template <class data_t>
template <template <typename> class vars_t>
auto Cell<data_t>::load_vars() const
{
    vars_t<data_t> vars;
    load_vars(vars);
    return vars;
}

template <class data_t>
ALWAYS_INLINE void Cell<data_t>::store_vars(const data_t &value,
                                            const int ivar) const
{
    SIMDIFY<data_t>(m_box_pointers.m_out_ptr[ivar])[m_out_index] = value;
}

template <class data_t>
template <int start_var, int end_var>
ALWAYS_INLINE void Cell<data_t>::store_vars(
    const tensor<1, data_t, GRInterval<start_var, end_var>::size()> &values,
    GRInterval<start_var, end_var> interval) const
{
    for (int i = 0; i < interval.size(); ++i)
        store_vars(values[i], interval.begin() + i);
}

template <class data_t>
ALWAYS_INLINE void
Cell<data_t>::store_vars(const std::array<data_t, NUM_VARS> &values) const
{
    for (int ivar = 0; ivar < NUM_VARS; ++ivar)
    {
        store_vars(values[ivar], ivar);
    }
}

template <class data_t>
template <template <typename> class vars_t>
void Cell<data_t>::store_vars(vars_t<data_t> &vars) const
{
    vars.enum_mapping([&](const int &ivar, data_t &var) {
        SIMDIFY<data_t>(m_box_pointers.m_out_ptr[ivar])[m_out_index] = var;
    });
}

#endif /* CELL_IMPL_HPP_ */
