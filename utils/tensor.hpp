#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "GRUtils.H"
#include "always_inline.hpp"
#include <initializer_list>

// C++ standard, 12.8 Copying class objects:
// Each subobject is assigned in the manner appropriate to its type:
// - if the subobject is an array, EACH ELEMENT IS ASSIGNED, in the manner appropriate to the element type;
template <int N, class data_t>
class tensor
{
    template <int, class> friend class tensor;
    typedef typename tensor<N-1, data_t>::arr_t arr_t[IDX_SPACEDIM];
    arr_t arr;

public:
    ALWAYS_INLINE
    tensor()
    {}

    ALWAYS_INLINE
    tensor(std::initializer_list<data_t> list) :
        arr (list)
    {}

    operator arr_t& ()
    {
        return arr;
    }

    operator const arr_t& () const
    {
        return arr;
    }
};

template <class data_t>
class tensor<0, data_t>
{
    template <int, class> friend class tensor;
    typedef data_t arr_t;
    arr_t arr;

public:
    ALWAYS_INLINE
    tensor()
    {}
    
    operator arr_t& ()
    {
        return arr;
    }

    operator const arr_t& () const
    {
        return arr;
    }
};

#endif /* TENSOR_HPP_ */
