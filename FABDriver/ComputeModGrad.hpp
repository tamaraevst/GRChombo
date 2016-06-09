//This class computes the modulus of the magnitude of chi
#ifndef COMPUTEMODGRAD_HPP_
#define COMPUTEMODGRAD_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"

#include <array>

class ComputeModGrad
{
protected:
   const FABDriverBase& m_driver;
   const FourthOrderDerivatives m_deriv;

public:
   ComputeModGrad(const FABDriverBase& driver, double dx) :
      m_driver (driver),
      m_deriv (dx, m_driver)
   {};

   template <class data_t>
      void compute(int x, int y, int z)
      {
         idx_t<data_t> idx = m_driver.in_idx(x, y, z);

         std::array<data_t, c_NUM> d1_arr[CH_SPACEDIM]; //Derivative index first
         for (int dir = 0; dir < CH_SPACEDIM; ++dir)
         {
            d1_arr[dir] = m_deriv.diff1(idx, dir);
         }

         std::array<data_t, c_NUM> mod_d1_arr = {};
         for (int comp = 0; comp < c_NUM; ++comp)
         {
            for (int dir = 0; dir < CH_SPACEDIM; ++dir)
            {
               mod_d1_arr[comp] += d1_arr[dir][comp]*d1_arr[dir][comp];
            }
            mod_d1_arr[comp] = simd_sqrt(mod_d1_arr[comp]);
         }

         // TODO: I really do not like this, but cannot think of a better way to do it yet...
         idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
         for (int comp = 0; comp < c_NUM; ++comp)
         {
            SIMDIFY<data_t>(m_driver.m_out_ptr[comp])[out_idx] = mod_d1_arr[comp];
         }
      }
   
};

#endif /* COMPUTEMODGRAD_HPP_ */
