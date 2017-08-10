#if !defined(CELL_HPP_)
#error "This file should only be included through Cell.hpp"
#endif

#ifndef CELL_IMPL_HPP_
#define CELL_IMPL_HPP_

#include "simd.hpp"
#include "Interval.H"

template <class data_t>
ALWAYS_INLINE
data_t
Cell<data_t>::local_vars(const int icomp) const
{
    return SIMDIFY<data_t>(m_box_pointers.m_in_ptr[icomp])[m_in_index];
}

template <class data_t>
ALWAYS_INLINE
void
Cell<data_t>::local_vars(data_t& out, const int icomp) const
{
    out = local_vars(icomp);
}

template <class data_t>
void
Cell<data_t>::local_vars(data_t (&out)[c_NUM]) const
{
    FORVARS(i) out[i] = local_vars(i);
}

template <class data_t>
template<template<typename> class vars_t>
void
Cell<data_t>::local_vars(vars_t<data_t>& vars) const
{
    vars.enum_mapping([&](const int& ivar, double& var)
                      { var = SIMDIFY<data_t>(m_box_pointers.m_in_ptr[ivar])[m_in_index]; });
}

template <class data_t>
template<template<typename> class vars_t>
vars_t<data_t>
Cell<data_t>::local_vars() const
{
    vars_t<data_t> vars;
    local_vars(vars);
    return vars;
}

template <class data_t>
ALWAYS_INLINE
void
Cell<data_t>::store_vars(const data_t& value, const int icomp) const
{
    SIMDIFY<data_t>(m_box_pointers.m_out_ptr[icomp])[m_out_index] = value;
}

template <class data_t>
template <int num_comp>
ALWAYS_INLINE
void
Cell<data_t>::store_vars(const tensor<1, data_t, num_comp>& values, const int start_comp) const
{
    for (int i = 0; i < num_comp; ++i) store_vars(values[i], start_comp + i);
}

template <class data_t>
ALWAYS_INLINE
void
Cell<data_t>::store_vars(const std::array<data_t, c_NUM>& values) const
{
    FORVARS(i) store_vars(values[i], i);
}

///This function stores all variables that have a corresponding value in a VarsBase object.
template <class data_t>
template<template<typename> class vars_t>
void
Cell<data_t>::store_vars(vars_t<data_t>& vars) const
{
    vars.enum_mapping([&](const int& ivar, double& var)
                      { SIMDIFY<data_t>(m_box_pointers.m_out_ptr[ivar])[m_out_index] = var; });
}

#endif /* CELL_IMPL_HPP_ */
