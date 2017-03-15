#if !defined(FABDRIVERBASE_HPP_)
#error "This file should only be included through FABDriverBase.hpp"
#endif

#ifndef FABDRIVERBASE_IMPL_HPP_
#define FABDRIVERBASE_IMPL_HPP_

#include "simd.hpp"
#include "Interval.H"
#include "CellIndex.hpp"

//Note: the compiler will not be able to deduce the template type.
//The usage therefore has to be local_vars<whatever data type you want>
template <class data_t>
ALWAYS_INLINE
data_t
FABDriverBase::local_vars(CellIndexIn in_index, const int icomp) const
{
    return SIMDIFY<data_t>(m_in_ptr[icomp])[in_index];
}

//This has the same functionality as above but allows the compiler to
//do type deduction on for data_t on out ... depends on taste.
template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::local_vars(data_t& out, CellIndexIn in_index, const int icomp) const
{
    out = SIMDIFY<data_t>(m_in_ptr[icomp])[in_index];
}

template <class data_t>
void
FABDriverBase::local_vars(data_t (&out)[c_NUM], CellIndexIn in_index) const
{
    FORVARS(i) out[i] = local_vars<data_t>(i, in_index);
}

template <class data_t>
void
FABDriverBase::local_vars(VarsBase<data_t>& vars, CellIndexIn in_index) const
{
    FORVARS(i) vars.assign(SIMDIFY<data_t>(m_in_ptr[i])[in_index], i);
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const data_t& value, CellIndexOut out_index, const int icomp) const
{
    SIMDIFY<data_t>(m_out_ptr[icomp])[out_index] = value;
}

template <class data_t>
ALWAYS_INLINE
void
FABDriverBase::store_vars(const std::array<data_t, c_NUM>& values, CellIndexOut out_index) const
{
    FORVARS(i) store_vars(values[i], out_index, i);
}

template <class data_t>
void
FABDriverBase::store_vars(const VarsBase<data_t>& vars, CellIndexOut out_index, const Interval a_comps) const
{
    for (int icomp=a_comps.begin(); icomp<=a_comps.end(); ++icomp)
    {
        CH_assert(vars.variable_defined(icomp));
        store_vars<data_t>(vars.template read<data_t>(icomp), out_index, icomp);
    }
}

template <class data_t>
void
//This function stores all variables that have a corresponding value in a VarsBase object.
//It will cycle through all components and check whether they have been assigned in the vars object.
//Avoid use for vars objects that only contain few components.
FABDriverBase::store_vars(const VarsBase<data_t>& vars, CellIndexOut out_index) const
{
    FORVARS(i)
    {
        //Only store the variable if it is defined in vars. Otherwise don't touch.
        if ( vars.variable_defined(i) ) store_vars<data_t>(vars.template read<data_t>(i), out_index, i);
    }
}

#endif /* FABDRIVERBASE_IMPL_HPP_ */
