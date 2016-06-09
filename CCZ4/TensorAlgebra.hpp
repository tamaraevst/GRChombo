#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "always_inline.hpp"
#include "tensor.hpp"

template <class data_t>
ALWAYS_INLINE
data_t
compute_trace(const tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &inverse_metric)
{
   data_t trace = 0;
   FOR2(i,j)
   {
      trace += inverse_metric[i][j]*tensor_LL[i][j];
   }
   return trace;
}

template <class data_t>
ALWAYS_INLINE
void
make_trace_free(tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &metric, const tensor<2,data_t> &inverse_metric)
{
   auto trace = compute_trace(tensor_LL, inverse_metric);
   FOR2(i,j)
   {
      tensor_LL[i][j] += - 1./(GR_SPACEDIM-1.) * metric[i][j] * trace;
   }
}

template <class data_t>
ALWAYS_INLINE
tensor<2,data_t>
raise_all(const tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &inverse_metric)
{
   tensor<2, data_t> tensor_UU;
   FOR2(i,j)
   {
      tensor_UU[i][j] = 0;
      FOR2(k,l)
      {
         tensor_UU[i][j] += inverse_metric[i][k]*inverse_metric[j][l]*tensor_LL[k][l];
      }
   }
   return tensor_UU;
}

#endif /* TENSORALGEBRA_HPP_ */
