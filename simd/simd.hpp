#ifndef SIMD_HPP_
#define SIMD_HPP_

#include <cmath>
#include "always_inline.hpp"

template <typename t> struct _simd_remove_const { typedef t type; };
template <typename t> struct _simd_remove_const<const t> { typedef t type; };

template <typename t> struct _simd_remove_volatile { typedef t type; };
template <typename t> struct _simd_remove_volatile<volatile t> { typedef t type; };

template <typename t> struct _simd_remove_cv {
    typedef typename _simd_remove_const<typename _simd_remove_volatile<t>::type>::type type;
};

template <typename q1, typename q2, typename t> struct _simd_enable_if_same { };
template <typename q, typename t> struct _simd_enable_if_same<q, q, t> { typedef t type; };

//This struct can be used to switch between simd and non-simd versions of the same function by overloading
struct disable_simd {};

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
    
    ALWAYS_INLINE
    t operator[] (int index) const { return m_value; }
    
    template <typename op_t>
    ALWAYS_INLINE
    simd foreach(op_t op) const
    {
        return simd(op(m_value));
    }
    
    template <typename op_t>
    ALWAYS_INLINE
    simd foreach(op_t op, t arg) const
    {
        return simd(op(m_value,arg));
    }
};

#define define_simd_overload(op) \
template <typename t> \
ALWAYS_INLINE \
simd<t> op(const simd<t>& a)\
{\
   return a.foreach(([&](t x) { return op(x);}));\
}

#define define_binary_simd_overload(op) \
template <typename t> \
ALWAYS_INLINE \
simd<t> op(const simd<t>& a, const simd<t>&b)\
{\
   return a.foreach(([&](t x, t arg) { return op(x,arg);}),b);\
}

/* Trascendental support:                               */
/* exp, sin, cos, log, sqrt, pow, tanh, tan, sinh, cosh */

define_simd_overload(exp)
define_simd_overload(exp2)
define_simd_overload(sin)
define_simd_overload(cos)
define_simd_overload(log)
define_simd_overload(log2)
define_simd_overload(sqrt)
define_binary_simd_overload(pow)
define_simd_overload(abs)
define_simd_overload(tanh)
define_simd_overload(tan)
define_simd_overload(sinh)
define_simd_overload(cosh)
define_simd_overload(acos)
define_simd_overload(asin)
define_simd_overload(atan)
define_binary_simd_overload(atan2)

/* Extra pow overloads */
template <typename t, typename t1> 
ALWAYS_INLINE 
simd<t> pow(const simd<t>& a, const t1 b)\
{
   simd<t> simd_b(b);
   return pow(a,simd_b);\
}

/* Extra atan2 overloads */
template <typename t, typename t1> 
ALWAYS_INLINE 
simd<t> atan2(const t1 b, const simd<t>& a)\
{
   simd<t> simd_b(b);
   return atan2(simd_b, a);\
}

#include "simd_base.hpp"

#if defined(__x86_64__)
#include "x64/x64.hpp"
#endif

template <typename t>
ALWAYS_INLINE
ostream& operator<< (ostream& os, const simd<t>& in_simd)
{
    t in_arr[simd_traits<t>::simd_len];
    simd<t>::store(in_arr, in_simd);

    os << "( ";
    for (int i = 0; i < simd_traits<t>::simd_len; ++i)
    {
        os << in_simd[i] << " ";
    }
    os << ")";
    if (os.fail())
        MayDay::Error("operator<<(ostream&,simd<t>&) failed");
    return os;
}

#include "simdify.hpp"

#endif /* SIMD_HPP_ */
