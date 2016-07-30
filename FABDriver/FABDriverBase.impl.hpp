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
std::array<data_t, c_NUM>
FABDriverBase::local_vars(idx_t<data_t> idx) const
{
    std::array<data_t, c_NUM> out;
    for (int i = 0; i < c_NUM; ++i)
    {
        out[i] = SIMDIFY<data_t>(m_in_ptr[i])[idx];
    }
    return out;
}

template <class data_t>
void
FABDriverBase::local_vars(VarsBase<data_t>& vars, idx_t<data_t> idx) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        vars.assign(SIMDIFY<data_t>(m_in_ptr[i])[idx], i);
    }
}

#endif /* FABDRIVERBASE_IMPL_HPP_ */
