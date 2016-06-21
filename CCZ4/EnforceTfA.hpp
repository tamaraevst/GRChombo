//This class enforces A to be trace-free
#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"

#include <array>

class EnforceTfA
{
protected:
   const FABDriverBase& m_driver;

public:
   EnforceTfA(const FABDriverBase& driver) :
      m_driver (driver)
   {}

   template <class data_t>
   struct vars_t
   {
      tensor<2, data_t> h;
      tensor<2, data_t> A;

      vars_t(){}

      template <class arr_t>
      vars_t(const arr_t& in)
      {
         h[0][0]  = in[c_h11];
         h[0][1]  = in[c_h12];
         h[0][2]  = in[c_h13];
         h[1][1]  = in[c_h22];
         h[1][2]  = in[c_h23];
         h[2][2]  = in[c_h33];

         h[1][0] = h[0][1];
         h[2][0] = h[0][2];
         h[2][1] = h[1][2];

         A[0][0]  = in[c_A11];
         A[0][1]  = in[c_A12];
         A[0][2]  = in[c_A13];
         A[1][1]  = in[c_A22];
         A[1][2]  = in[c_A23];
         A[2][2]  = in[c_A33];

         A[1][0] = A[0][1];
         A[2][0] = A[0][2];
         A[2][1] = A[1][2];
      }
   };

   template <class data_t>
   void compute(int x, int y, int z)
   {
      idx_t<data_t> idx = m_driver.in_idx(x, y, z);

      vars_t<data_t> vars = m_driver.local_vars(idx);

      auto h_UU = TensorAlgebra::compute_inverse(vars.h);
      TensorAlgebra::make_trace_free(vars.A, vars.h, h_UU);

      idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A11])[out_idx] = vars.A[0][0];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A12])[out_idx] = vars.A[0][1];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A13])[out_idx] = vars.A[0][2];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A22])[out_idx] = vars.A[1][1];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A23])[out_idx] = vars.A[1][2];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_A33])[out_idx] = vars.A[2][2];
   }
};

#endif /* FIXTFA_HPP_ */
