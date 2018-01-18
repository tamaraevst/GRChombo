#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"

/// This class implements a tensor with given rank, element data type, and
/// dimension.  By default the dimension is equal to DEFAULT_TENSOR_DIM.
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
