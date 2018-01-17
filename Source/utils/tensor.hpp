#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "GRUtils.hpp"
#include "AlwaysInline.hpp"
//#include <initializer_list>

// C++ standard, 12.8 Copying class objects:
// Each subobject is assigned in the manner appropriate to its type:
// - if the subobject is an array, EACH ELEMENT IS ASSIGNED, in the manner
// appropriate to the element type;
template <int rank, class data_t, int size = DEFAULT_TENSOR_DIM> class tensor
{
    template <int, class, int> friend class tensor;
    typedef typename tensor<rank - 1, data_t>::arr_t arr_t[size];
    arr_t arr;

  public:
    ALWAYS_INLINE
    tensor() {}

    //    ALWAYS_INLINE
    //    tensor(std::initializer_list<data_t> list) :
    //        arr (list)
    //    {}

    template <typename... T> ALWAYS_INLINE tensor(T... data) : arr{data...} {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

template <class data_t, int size> class tensor<0, data_t, size>
{
    template <int, class, int> friend class tensor;
    typedef data_t arr_t;
    arr_t arr;

  public:
    ALWAYS_INLINE
    tensor() {}

    ALWAYS_INLINE
    tensor(data_t val) : arr(val) {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

#endif /* TENSOR_HPP_ */
