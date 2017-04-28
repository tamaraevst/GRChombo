#if !defined(CELL_HPP_)
#error "This file should only be included through Cell.hpp"
#endif

#ifndef CELL_IMPL_HPP_
#define CELL_IMPL_HPP_

#include "simd.hpp"
#include "Interval.H"

//Note: the compiler will not be able to deduce the template type.
//The usage therefore has to be local_vars<whatever data type you want>
template <class data_t>
ALWAYS_INLINE
data_t
Cell::local_vars(const int icomp) const
{
    return SIMDIFY<data_t>(m_box_pointers.m_in_ptr[icomp])[m_in_index];
}

//This has the same functionality as above but allows the compiler to
//do type deduction on for data_t on out ... depends on taste.
template <class data_t>
ALWAYS_INLINE
void
Cell::local_vars(data_t& out, const int icomp) const
{
    out = SIMDIFY<data_t>(m_box_pointers.m_in_ptr[icomp])[m_in_index];
}

template <class data_t>
void
Cell::local_vars(data_t (&out)[c_NUM]) const
{
    FORVARS(i) out[i] = local_vars<data_t>(i);
}

template <class data_t>
void
Cell::local_vars(VarsBase<data_t>& vars) const
{
    FORVARS(i) vars.assign(SIMDIFY<data_t>(m_box_pointers.m_in_ptr[i])[m_in_index], i);
}

template <class data_t>
ALWAYS_INLINE
void
Cell::store_vars(const data_t& value, const int icomp) const
{
    SIMDIFY<data_t>(m_box_pointers.m_out_ptr[icomp])[m_out_index] = value;
}

template <class data_t>
ALWAYS_INLINE
void
Cell::store_vars(const std::array<data_t, c_NUM>& values) const
{
    FORVARS(i) store_vars(values[i], i);
}

template <class data_t>
void
Cell::store_vars(const VarsBase<data_t>& vars, const Interval a_comps) const
{
    for (int icomp=a_comps.begin(); icomp<=a_comps.end(); ++icomp)
    {
        CH_assert(vars.variable_defined(icomp));
        store_vars<data_t>(vars.template read<data_t>(icomp), icomp);
    }
}

template <class data_t>
void
//This function stores all variables that have a corresponding value in a VarsBase object.
//It will cycle through all components and check whether they have been assigned in the vars object.
//Avoid use for vars objects that only contain few components.
Cell::store_vars(const VarsBase<data_t>& vars) const
{
    FORVARS(i)
    {
        //Only store the variable if it is defined in vars. Otherwise don't touch.
        if ( vars.variable_defined(i) ) store_vars<data_t>(vars.template read<data_t>(i), i);
    }
}

#endif /* CELL_IMPL_HPP_ */
