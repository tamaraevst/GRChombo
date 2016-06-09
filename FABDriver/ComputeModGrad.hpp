//This class computes the modulus of the magnitude of chi
#ifndef COMPUTEMODGRAD_HPP_
#define COMPUTEMODGRAD_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"

class ComputeModGrad
{
public:
   ComputeModGrad(const FABDriverBase& driver, double dx) : m_dx (dx), m_driver(driver){};

   template <class data_t>
      void compute(int x, int y, int z)
      {
         const int idx = m_driver.m_stride[2]*(z-m_driver.m_in_lo[2]) + m_driver.m_stride[1]*(y-m_driver.m_in_lo[1]) + (x-m_driver.m_in_lo[0]);

         data_t d1Arr[CH_SPACEDIM][c_NUM]; //Derivative index first
         for (int dir = 0; dir < CH_SPACEDIM; ++dir)
         {
            m_driver.diff1(idx, m_driver.m_stride[dir], m_dx, d1Arr[dir]);
         }

         data_t modD1Arr[c_NUM];
         for (int comp = 0; comp < c_NUM; ++comp)
         {
            modD1Arr[comp] = 0;
            for (int dir = 0; dir < CH_SPACEDIM; ++dir)
            {
               modD1Arr[comp] += d1Arr[dir][comp]*d1Arr[dir][comp];
            }
            modD1Arr[comp] = simd_sqrt(modD1Arr[comp]);
         }

         // TODO: I really do not like this, but cannot think of a better way to do it yet...
         const int out_idx = m_driver.m_out_stride[2]*(z-m_driver.m_out_lo[2]) + m_driver.m_out_stride[1]*(y-m_driver.m_out_lo[1]) + (x-m_driver.m_out_lo[0]);
         for (int comp = 0; comp < c_NUM; ++comp)
         {
            SIMDIFY<data_t>(m_driver.m_out_ptr[comp])[out_idx]    =  modD1Arr[comp];
         }
      }

protected:
   const double m_dx;
   const FABDriverBase& m_driver;
};

#endif /* COMPUTEMODGRAD_HPP_ */
