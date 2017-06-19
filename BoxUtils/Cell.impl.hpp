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
void
Cell<data_t>::local_vars(VarsBase<data_t>& vars) const
{
    FORVARS(i) vars.assign(SIMDIFY<data_t>(m_box_pointers.m_in_ptr[i])[m_in_index], i);
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

template <class data_t>
void
Cell<data_t>::store_vars(const VarsBase<data_t>& vars, const Interval a_comps) const
{
    for (int icomp=a_comps.begin(); icomp<=a_comps.end(); ++icomp)
    {
        CH_assert(vars.variable_defined(icomp));
        store_vars(vars.template read<data_t>(icomp), icomp);
    }
}

///This function stores all variables that have a corresponding value in a VarsBase object.
/**It will cycle through all components and check whether they have been assigned in the vars object.
  *Avoid use for vars objects that only contain few components.
  */
template <class data_t>
void
Cell<data_t>::store_vars(const VarsBase<data_t>& vars) const
{
    FORVARS(i)
    {
        //Only store the variable if it is defined in vars. Otherwise don't touch.
        if ( vars.variable_defined(i) ) store_vars(vars.template read<data_t>(i), i);
    }
}

#endif /* CELL_IMPL_HPP_ */
