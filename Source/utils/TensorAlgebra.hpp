#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "tensor.hpp"

template <class data_t> struct chris_t
{
    tensor<3, data_t> ULL;        //!<standard christoffel symbols
    tensor<3, data_t> LLL;        //!<3 lower indices
    tensor<1, data_t> contracted; //!<contracted christoffel
};

namespace TensorAlgebra
{
/// Computes determinant of a symmetric 3x3 matrix
template <class data_t>
ALWAYS_INLINE data_t compute_determinant_sym(const tensor<2, data_t, 3> &matrix)
{
    data_t det = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                 2 * matrix[0][1] * matrix[0][2] * matrix[1][2] -
                 matrix[0][0] * matrix[1][2] * matrix[1][2] -
                 matrix[1][1] * matrix[0][2] * matrix[0][2] -
                 matrix[2][2] * matrix[0][1] * matrix[0][1];

    return det;
}

/// Computes the determinant of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
ALWAYS_INLINE data_t compute_determinant(const tensor<2, data_t, 3> &matrix)
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

/// Computes the inverse of a symmetric 3x3 matrix directly.
template <class data_t>
tensor<2, data_t> compute_inverse_sym(const tensor<2, data_t, 3> &matrix)
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

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
template <class data_t>
tensor<2, data_t> compute_inverse(const tensor<2, data_t, 3> &matrix)
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

/// Computes the trace of a 2-tensor with lower indices given an inverse metric.
template <class data_t>
ALWAYS_INLINE data_t compute_trace(const tensor<2, data_t> &tensor_LL,
                                   const tensor<2, data_t> &inverse_metric)
{
    data_t trace = 0.;
    FOR2(i, j) { trace += inverse_metric[i][j] * tensor_LL[i][j]; }
    return trace;
}

/// Computes the trace of a 1,1 tensor (a matrix) - no metric required.
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

/// Computes dot product of a vector and a covector (no metric required)
template <class data_t>
ALWAYS_INLINE data_t compute_dot_product(const tensor<1, data_t> &vector_U,
                                         const tensor<1, data_t> &covector_L)
{
    data_t dot_product = 0.;
    FOR1(i) dot_product += vector_U[i] * covector_L[i];
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
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

/// Removes the trace of a 2-tensor with lower indices given a metric and an
/// inverse metric.  Or a tensor with upper indices given an inverse metric and
/// a metric.
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

/// Raises the index of a covector
template <class data_t>
ALWAYS_INLINE tensor<1, data_t>
raise_all(const tensor<1, data_t> &tensor_L,
          const tensor<2, data_t> &inverse_metric)
{
    tensor<1, data_t> tensor_U = 0.;
    FOR2(i, j) { tensor_U[i] += inverse_metric[i][j] * tensor_L[j]; }
    return tensor_U;
}

/// Raises the indices of a 2-tensor
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

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readibility
template <class data_t>
ALWAYS_INLINE tensor<1, data_t> lower_all(const tensor<1, data_t> &tensor_U,
                                          const tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-tensor
/// Note: same functionality as raise; included to improve readibility
template <class data_t>
ALWAYS_INLINE tensor<2, data_t> lower_all(const tensor<2, data_t> &tensor_UU,
                                          const tensor<2, data_t> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Computes the (i,j) component of the Kronecker delta
constexpr int delta(int i, int j) { return (i == j); }

/// Computes the levi-civita symbol (3D, NB, symbol, not the tensor)
inline tensor<3,double> epsilon()
{
    tensor<3, double> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    return epsilon;
}

/// Computes the conformal christoffel symbol
template <class data_t>
chris_t<data_t>
compute_christoffel(const tensor<2, tensor<1, data_t>> &d1_metric,
                    const tensor<2, data_t> &h_UU)
{
    chris_t<data_t> out;

    FOR3(i, j, k)
    {
        out.LLL[i][j][k] = 0.5 * (d1_metric[j][i][k] + d1_metric[k][i][j] -
                                  d1_metric[j][k][i]);
    }
    FOR3(i, j, k)
    {
        out.ULL[i][j][k] = 0;
        FOR1(l) { out.ULL[i][j][k] += h_UU[i][l] * out.LLL[l][j][k]; }
    }
    FOR1(i)
    {
        out.contracted[i] = 0;
        FOR2(j, k) { out.contracted[i] += h_UU[j][k] * out.ULL[i][j][k]; }
    }

    return out;
}
}

#endif /* TENSORALGEBRA_HPP_ */
