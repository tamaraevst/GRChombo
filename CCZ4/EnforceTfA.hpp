//This class enforces A to be trace-free
#ifndef FIXTFA_HPP_
#define FIXTFA_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "CCZ4Geometry.hpp"
#include "TensorAlgebra.hpp"

class EnforceTfA
{
   public:
      EnforceTfA(const FABDriverBase& driver) : m_driver(driver){};

      template <class data_t>
      struct vars_t
      {
         tensor<2, data_t> h;
         tensor<2, data_t> A;
      };

      template <class data_t>
      void compute(int x, int y, int z)
      {
         const int idx = m_driver.m_stride[2]*(z-m_driver.m_in_lo[2]) + m_driver.m_stride[1]*(y-m_driver.m_in_lo[1]) + (x-m_driver.m_in_lo[0]);

         vars_t<data_t> vars;
         {
            data_t varsArr[c_NUM];
            m_driver.local_vars(idx, varsArr);
            demarshall(varsArr, vars);
         }
         auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);
         make_trace_free(vars.A, vars.h, h_UU);
      }

   protected:
      const FABDriverBase& m_driver;

      template <class data_t>
      void demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out)
      {
         out.h[0][0]  = in[c_h11];
         out.h[0][1]  = in[c_h12];
         out.h[0][2]  = in[c_h13];
         out.h[1][1]  = in[c_h22];
         out.h[1][2]  = in[c_h23];
         out.h[2][2]  = in[c_h33];

         out.h[1][0] = out.h[0][1];
         out.h[2][0] = out.h[0][2];
         out.h[2][1] = out.h[1][2];

         out.A[0][0]  = in[c_A11];
         out.A[0][1]  = in[c_A12];
         out.A[0][2]  = in[c_A13];
         out.A[1][1]  = in[c_A22];
         out.A[1][2]  = in[c_A23];
         out.A[2][2]  = in[c_A33];

         out.A[1][0] = out.A[0][1];
         out.A[2][0] = out.A[0][2];
         out.A[2][1] = out.A[1][2];
      }
};

#endif /* FIXTFA_HPP_ */
