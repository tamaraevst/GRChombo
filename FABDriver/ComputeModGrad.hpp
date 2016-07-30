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

       tensor<1,data_t> d1_arr[c_NUM];
       FOR1(idir) m_deriv.diff1(d1_arr,idx, idir);

       std::array<data_t, c_NUM> mod_d1_arr = {0.};
       FORVARS(ivar)
       {
           FOR1(idir)
           {
               mod_d1_arr[ivar] += d1_arr[ivar][idir]*d1_arr[ivar][idir];
           }
           mod_d1_arr[ivar] = simd_sqrt(mod_d1_arr[ivar]);
       }

       // TODO: I really do not like this, but cannot think of a better way to do it yet...
       idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
       FORVARS(ivar)
       {
           SIMDIFY<data_t>(m_driver.m_out_ptr[ivar])[out_idx] = mod_d1_arr[ivar];
       }
   }

};

#endif /* COMPUTEMODGRAD_HPP_ */
