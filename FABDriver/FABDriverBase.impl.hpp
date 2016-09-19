#if !defined(FABDRIVERBASE_HPP_)
#error "This file should only be included through FABDriverBase.hpp"
#endif

#ifndef FABDRIVERBASE_IMPL_HPP_
#define FABDRIVERBASE_IMPL_HPP_

#include "simd.hpp"

ALWAYS_INLINE
int
FABDriverBase::in_idx(int ix, int iy, int iz) const
{
    return m_in_stride[2]*(iz-m_in_lo[2]) + m_in_stride[1]*(iy-m_in_lo[1]) + (ix-m_in_lo[0]);
}

ALWAYS_INLINE
int
FABDriverBase::out_idx(int ix, int iy, int iz) const
{
    return m_out_stride[2]*(iz-m_out_lo[2]) + m_out_stride[1]*(iy-m_out_lo[1]) + (ix-m_out_lo[0]);
}

template <class data_t>
ALWAYS_INLINE
data_t
FABDriverBase::local_vars(const idx_t<data_t> idx, const int icomp) const
{
    return SIMDIFY<data_t>(m_in_ptr[icomp])[idx];
}

template <class data_t>
std::array<data_t, c_NUM>
FABDriverBase::local_vars(const idx_t<data_t> idx) const
{
    std::array<data_t, c_NUM> out;
    FORVARS(i) out[i] = local_vars(idx, i);
    return out;
}

template <class data_t>
void
FABDriverBase::local_vars(VarsBase<data_t>& vars, const idx_t<data_t> idx) const
{
    FORVARS(i) vars.assign(SIMDIFY<data_t>(m_in_ptr[i])[idx], i);
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const data_t& value, const idx_t<data_t> out_idx, const int icomp) const
{
    SIMDIFY<data_t>(m_out_ptr[icomp])[out_idx] = value;
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const std::array<data_t, c_NUM>& values, const idx_t<data_t> out_idx) const
{
    FORVARS(i) store_vars(values[i], out_idx, i);
}

template <class data_t>
void
//This function stores all variables that have a corresponding value in a VarsBase object.
//It will cycle through all components and check whether they have been assigned in the vars object.
//Avoid use for vars objects that only contain few components.
FABDriverBase::store_vars(const VarsBase<data_t>& vars, const idx_t<data_t> out_idx) const
{
    FORVARS(i)
    {
        //Only store the variable if it is defined in vars. Otherwise don't touch.
        if ( vars.variable_defined(i) ) store_vars(vars.template read<data_t>(i), out_idx, i);
    }
}

#endif /* FABDRIVERBASE_IMPL_HPP_ */
