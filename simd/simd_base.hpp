#ifndef SIMD_BASE_HPP_
#define SIMD_BASE_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

#include <cmath>
#include <type_traits>

template <typename t>
struct simd_traits
{
    typedef t data_t;
    typedef bool mask_t;
    static const int simd_len = 1;
};

template <typename t>
struct simd_base
{
    static const int simd_len = simd_traits<t>::simd_len;
    typedef typename simd_traits<t>::data_t simd_data_t;
    simd_data_t m_value;

    ALWAYS_INLINE
    simd_base() :
        m_value ()
    {}

#pragma GCC diagnostic ignored "-Wuninitialized"
    ALWAYS_INLINE
    simd_base(const simd_data_t& x) :
        m_value (x)
    {}
#pragma GCC diagnostic pop

    ALWAYS_INLINE
    operator simd_data_t&()
    {
        return m_value;
    }

    ALWAYS_INLINE
    operator const simd_data_t&() const
    {
        return m_value;
    }

    // Note: These binary ops allow us to write e.g. simd<double>+int
    // or 2*y etc. The "friend function" construction makes
    // sure that these NON-MEMBER binary operations
    // are instantiated AT THE SAME TIME as simd<t>.
    // The alternative is to template
    // on the two operand types, however this gives rise to
    // headache-inducing ambiguity because we have defined
    // casts between t and simd<t>
    friend ALWAYS_INLINE
    simd<t> operator+(const simd<t>& a, const simd<t>& b)
    {
        simd<t> out(static_cast<const simd<t>&>(a));
        out += static_cast<const simd<t>&>(b);
        return out;
    }

    friend ALWAYS_INLINE
    simd<t> operator-(const simd<t>& a, const simd<t>& b)
    {
        simd<t> out(static_cast<const simd<t>&>(a));
        out -= static_cast<const simd<t>&>(b);
        return out;
    }

    friend ALWAYS_INLINE
    simd<t> operator-(const simd<t>& a)
    {
        simd<t> out(0);
        out -= static_cast<const simd<t>&>(a);
        return out;
    }

    friend ALWAYS_INLINE
    simd<t> operator*(const simd<t>& a, const simd<t>& b)
    {
        simd<t> out(static_cast<const simd<t>&>(a));
        out *= static_cast<const simd<t>&>(b);
        return out;
    }

    friend ALWAYS_INLINE
    simd<t> operator/(const simd<t>& a, const simd<t>& b)
    {
        simd<t> out(static_cast<const simd<t>&>(a));
        out /= static_cast<const simd<t>&>(b);
        return out;
    }
};

template <typename t> struct simd_compare { typedef bool type; };
template <typename t> struct simd_compare<simd<t> > { typedef typename simd_traits<t>::mask_t type; };

template <typename t>
t simd_conditional(typename simd_compare<t>::type cond, const t& true_value, const t& false_value)
{
    return cond ? true_value : false_value;
}

template <typename t>
ALWAYS_INLINE
t simd_min(const t& a, const t& b)
{
    return (a <= b) ? a : b;
}

template <typename t>
ALWAYS_INLINE
t simd_max(const t& a, const t& b)
{
    return (a > b) ? a : b;
}

template <typename t>
ALWAYS_INLINE
typename simd_compare<t>::type
simd_compare_lt(const t& a, const t& b)
{
    return a < b;
}

template <typename t>
ALWAYS_INLINE
typename simd_compare<t>::type
simd_compare_gt(const t& a, const t& b)
{
    return a > b;
}

template <typename t>
ALWAYS_INLINE
t simd_sqrt(const t& a)
{
   return sqrt(a);
}

template <typename t>
ALWAYS_INLINE
t simd_log(const t& a)
{
   return log(a);
}

//This function will only be enabled if simd<t> is not the same as t, i.e. if we could potentially vectorise
template <typename t>
ALWAYS_INLINE
typename std::enable_if<(simd_traits<t>::simd_len > 1), simd<t> >::type
simd_log(const simd<t>& a)
{
   //For log we just let the compiler auto-vectorise
   t in_arr[simd_traits<t>::simd_len];
   t out_arr[simd_traits<t>::simd_len];
   simd<t>::store(in_arr, a);

   #pragma simd vectorlengthfor(t)
   for (int i = 0; i < simd_traits<t>::simd_len; ++i)
   {
      out_arr[i] = log(in_arr[i]);
   }

   return simd<t>::load(out_arr);
}

#endif /* SIMD_BASE_HPP_ */
