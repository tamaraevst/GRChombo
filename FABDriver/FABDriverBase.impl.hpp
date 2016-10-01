#if !defined(FABDRIVERBASE_HPP_)
#error "This file should only be included through FABDriverBase.hpp"
#endif

#ifndef FABDRIVERBASE_IMPL_HPP_
#define FABDRIVERBASE_IMPL_HPP_

#include "simd.hpp"
#include "Interval.H"

ALWAYS_INLINE
void
FABDriverBase::set_idx(int ix, int iy, int iz)
{
    m_in_idx  = m_in_stride[2]*(iz-m_in_lo[2]) + m_in_stride[1]*(iy-m_in_lo[1]) + (ix-m_in_lo[0]);
    m_out_idx = m_out_stride[2]*(iz-m_out_lo[2]) + m_out_stride[1]*(iy-m_out_lo[1]) + (ix-m_out_lo[0]);
}

ALWAYS_INLINE
int
FABDriverBase::get_in_idx() const
{
    return m_in_idx;
}

ALWAYS_INLINE
int
FABDriverBase::get_out_idx() const
{
    return m_out_idx;
}

//Note: the compiler will not be able to deduce the template type.
//The usage therefore has to be local_vars<whatever data type you want>
template <class data_t>
ALWAYS_INLINE
data_t
FABDriverBase::local_vars(const int icomp) const
{
    return SIMDIFY<data_t>(m_in_ptr[icomp])[m_in_idx];
}

//This has the same functionality as above but allows the compiler to
//do type deduction on for data_t on out ... depends on taste.
template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::local_vars(data_t& out, const int icomp) const
{
    out = SIMDIFY<data_t>(m_in_ptr[icomp])[m_in_idx];
}

template <class data_t>
void
FABDriverBase::local_vars(data_t (&out)[c_NUM]) const
{
    FORVARS(i) out[i] = local_vars<data_t>(i);
}

template <class data_t>
void
FABDriverBase::local_vars(VarsBase<data_t>& vars) const
{
    FORVARS(i) vars.assign(SIMDIFY<data_t>(m_in_ptr[i])[m_in_idx], i);
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const data_t& value, const int icomp) const
{
    SIMDIFY<data_t>(m_out_ptr[icomp])[m_out_idx] = value;
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const std::array<data_t, c_NUM>& values) const
{
    FORVARS(i) store_vars(values[i], i);
}

template <class data_t>
void
FABDriverBase::store_vars(const VarsBase<data_t>& vars, const Interval a_comps) const
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
FABDriverBase::store_vars(const VarsBase<data_t>& vars) const
{
    FORVARS(i)
    {
        //Only store the variable if it is defined in vars. Otherwise don't touch.
        if ( vars.variable_defined(i) ) store_vars<data_t>(vars.template read<data_t>(i), i);
    }
}

#endif /* FABDRIVERBASE_IMPL_HPP_ */
