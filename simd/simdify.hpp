#ifndef SIMDIFY_HPP_
#define SIMDIFY_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

template <typename data_t>
struct simd_proxy
{
    typedef simd<typename _simd_remove_cv<data_t>::type> simd_t;
    data_t *const m_ptr;

    ALWAYS_INLINE
    simd_proxy(data_t *m_ptr) :
        m_ptr (m_ptr)
    {}

    ALWAYS_INLINE
    data_t*
    operator&() const
    {
        return m_ptr;
    }
    
    ALWAYS_INLINE
    operator simd_t() const
    {
        return simd_t::load(m_ptr);
    }

    ALWAYS_INLINE
    simd_proxy& operator=(const simd_t& rhs)
    {
        simd_t::store(m_ptr, rhs);
        return *this;
    }

    /*ALWAYS_INLINE
    const simd_proxy& operator=(const simd_t& rhs) const
    {
        const_cast<simd_proxy&>(*this) = rhs;
    }*/

    ALWAYS_INLINE
    simd_t operator+(const simd_t& other)
    {
        return static_cast<simd_t>(*this) + other;
    }

    ALWAYS_INLINE
    simd_t operator-(const simd_t& other)
    {
        return static_cast<simd_t>(*this) - other;
    }

    ALWAYS_INLINE
    simd_t operator*(const simd_t& other)
    {
        return static_cast<simd_t>(*this) * other;
    }

    ALWAYS_INLINE
    simd_t operator/(const simd_t& other)
    {
        return static_cast<simd_t>(*this) / other;
    }

    friend ALWAYS_INLINE
    simd_t operator+(const simd_t& a, const simd_proxy& b)
    {
        return a + static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE
    simd_t operator-(const simd_t& a, const simd_proxy& b)
    {
        return a - static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE
    simd_t operator*(const simd_t& a, const simd_proxy& b)
    {
        return a * static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE
    simd_t operator/(const simd_t& a, const simd_proxy& b)
    {
        return a / static_cast<simd_t>(b);
    }

    ALWAYS_INLINE
    simd_proxy& operator+=(const simd_t& rhs)
    {
        (*this) = (*this) + rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy& operator-=(const simd_t& rhs)
    {
        (*this) = (*this) - rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy& operator*=(const simd_t& rhs)
    {
        (*this) = (*this) * rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy& operator/=(const simd_t& rhs)
    {
        (*this) = (*this) / rhs;
        return *this;
    }
};

template <typename data_t>
struct simd_array_wrapper
{
    data_t *m_ptr;

    ALWAYS_INLINE
    simd_array_wrapper(data_t *m_ptr) :
        m_ptr (m_ptr)
    {}

    ALWAYS_INLINE
    simd_proxy<data_t>
    operator[](int i) const
    {
        return &m_ptr[i];
    }

    ALWAYS_INLINE
    simd_proxy<data_t>
    operator*() const
    {
        return *this[0];
    }
};

template <typename t, typename ptr_t>
struct _simdify
{
    typedef typename _simd_enable_if_same<t, typename _simd_remove_cv<ptr_t>::type, ptr_t*>::type type;
};

template <typename t, typename ptr_t>
struct _simdify<simd<t>, ptr_t>
{
    typedef typename _simd_enable_if_same<t, typename _simd_remove_cv<ptr_t>::type, simd_array_wrapper<ptr_t> >::type type;
};

template <typename t, typename ptr_t>
ALWAYS_INLINE
typename _simdify<t, ptr_t>::type
SIMDIFY(ptr_t* ptr)
{
    return ptr;
}

#endif /* SIMDIFY_HPP_ */
