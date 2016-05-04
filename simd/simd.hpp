#ifndef SIMD_HPP_
#define SIMD_HPP_

#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((always_inline)) __inline__
#else
#define ALWAYS_INLINE inline
#endif

template <typename t> struct _simd_remove_const { typedef t type; };
template <typename t> struct _simd_remove_const<const t> { typedef t type; };

template <typename t> struct _simd_remove_volatile { typedef t type; };
template <typename t> struct _simd_remove_volatile<volatile t> { typedef t type; };

template <typename t> struct _simd_remove_cv {
    typedef typename _simd_remove_const<typename _simd_remove_volatile<t>::type>::type type;
};

template <typename q1, typename q2, typename t> struct _simd_enable_if_same { };
template <typename q, typename t> struct _simd_enable_if_same<q, q, t> { typedef t type; };

// Base template type: fallback for when there is no SIMD support for a data type
template <typename t>
struct simd
{
    t m_value;

    ALWAYS_INLINE
    simd() :
        m_value ()
    {}
    
    ALWAYS_INLINE
    simd(const t& x) :
        m_value (x)
    {}

    ALWAYS_INLINE
    operator t&()
    {
        return m_value;
    }

    ALWAYS_INLINE
    operator const t&() const
    {
        return m_value;
    }

    ALWAYS_INLINE
    static simd
    load(const double *ptr)
    {
        return *ptr;
    }

    ALWAYS_INLINE
    static void
    store(double *ptr, const simd& a)
    {
        *ptr = a.m_value;
    }

    ALWAYS_INLINE
    simd&
    operator+=(const simd& a)
    {
        m_value += a.m_value;
        return *this;
    }
    ALWAYS_INLINE
    simd&
    operator-=(const simd& a)
    {
        m_value -= a.m_value;
        return *this;
    }
    ALWAYS_INLINE
    simd&
    operator*=(const simd& a)
    {
        m_value *= a.m_value;
        return *this;
    }
    ALWAYS_INLINE
    simd&
    operator/=(const simd& a)
    {
        m_value /= a.m_value;
        return *this;
    }
};

#include "simd_base.hpp"

#if defined(__x86_64__)
#include "x64/x64.hpp"
#endif

#include "simdify.hpp"

#undef ALWAYS_INLINE

#endif /* SIMD_HPP_ */
