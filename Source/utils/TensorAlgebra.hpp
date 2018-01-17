#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "always_inline.hpp"
#include "tensor.hpp"

namespace TensorAlgebra
{
template <class data_t>
ALWAYS_INLINE data_t
compute_determinant_sym(const tensor<2, data_t, 3>
                            &matrix) // This function only works for 3D matrix
{
    data_t det = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                 2 * matrix[0][1] * matrix[0][2] * matrix[1][2] -
                 matrix[0][0] * matrix[1][2] * matrix[1][2] -
                 matrix[1][1] * matrix[0][2] * matrix[0][2] -
                 matrix[2][2] * matrix[0][1] * matrix[0][1];

    return det;
}

template <class data_t>
ALWAYS_INLINE data_t
compute_determinant(const tensor<2, data_t, 3>
                        &matrix) // This function only works for 3D matrix
{
    data_t det =
        matrix[0][0] *
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] *
            (matrix[2][2] * matrix[1][0] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] *
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return det;
}

template <class data_t>
tensor<2, data_t>
compute_inverse_sym(const tensor<2, data_t, 3>
                        &matrix) // This function only works for 3D matrix
{
    data_t deth = compute_determinant_sym(matrix);
    data_t deth_inverse = 1. / deth;
    tensor<2, data_t> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[1][2] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[0][1] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

template <class data_t>
tensor<2, data_t>
compute_inverse(const tensor<2, data_t, 3>
                    &matrix) // This function only works for 3D matrix
{
    data_t deth = compute_determinant(matrix);
    data_t deth_inverse = 1. / deth;
    tensor<2, data_t> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[2][0] * matrix[1][2] - matrix[1][0] * matrix[2][2]) *
                 deth_inverse;
    h_UU[1][0] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) *
                 deth_inverse;
    h_UU[2][0] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) *
                 deth_inverse;
    h_UU[2][1] = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;

    return h_UU;
}

template <class data_t>
ALWAYS_INLINE data_t compute_trace(const tensor<2, data_t> &tensor_LL,
                                   const tensor<2, data_t> &inverse_metric)
{
    data_t trace = 0.;
    FOR2(i, j) { trace += inverse_metric[i][j] * tensor_LL[i][j]; }
    return trace;
}

template <class data_t>
ALWAYS_INLINE data_t compute_trace(const tensor<2, data_t> &tensor_UL)
{
    data_t trace = 0.;
    FOR1(i) trace += tensor_UL[i][i];
    return trace;
}

template <class data_t>
ALWAYS_INLINE data_t
compute_trace(const tensor<1, tensor<1, data_t>> &tensor_UL)
{
    data_t trace = 0.;
    FOR1(i) trace += tensor_UL[i][i];
    return trace;
}

template <class data_t>
ALWAYS_INLINE data_t compute_dot_product(const tensor<1, data_t> &vector_U,
                                         const tensor<1, data_t> &covector_L)
{
    data_t dot_product = 0.;
    FOR1(i) dot_product += vector_U[i] * covector_L[i];
    return dot_product;
}

template <class data_t>
ALWAYS_INLINE data_t compute_dot_product(
    const tensor<1, data_t> &covector1_L, const tensor<1, data_t> &covector2_L,
    const tensor<2, data_t> &inverse_metric)
{
    data_t dot_product = 0.;
    FOR2(m, n)
    {
        dot_product += inverse_metric[m][n] * covector1_L[m] * covector2_L[n];
    }
    return dot_product;
}

template <class data_t>
ALWAYS_INLINE void make_trace_free(tensor<2, data_t> &tensor_LL,
                                   const tensor<2, data_t> &metric,
                                   const tensor<2, data_t> &inverse_metric)
{
    auto trace = compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR2(i, j)
    {
        tensor_LL[i][j] -= one_over_gr_spacedim * metric[i][j] * trace;
    }
}

template <class data_t>
ALWAYS_INLINE tensor<1, data_t>
raise_all(const tensor<1, data_t> &tensor_L,
          const tensor<2, data_t> &inverse_metric)
{
    tensor<1, data_t> tensor_U = 0.;
    FOR2(i, j) { tensor_U[i] += inverse_metric[i][j] * tensor_L[j]; }
    return tensor_U;
}

template <class data_t>
ALWAYS_INLINE tensor<2, data_t>
raise_all(const tensor<2, data_t> &tensor_LL,
          const tensor<2, data_t> &inverse_metric)
{
    tensor<2, data_t> tensor_UU = 0.;
    FOR4(i, j, k, l)
    {
        tensor_UU[i][j] +=
            inverse_metric[i][k] * inverse_metric[j][l] * tensor_LL[k][l];
    }
    return tensor_UU;
}

template <class data_t>
ALWAYS_INLINE tensor<1, data_t> lower_all(const tensor<1, data_t> &tensor_U,
                                          const tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

template <class data_t>
ALWAYS_INLINE tensor<2, data_t> lower_all(const tensor<2, data_t> &tensor_UU,
                                          const tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}
}

#endif /* TENSORALGEBRA_HPP_ */
