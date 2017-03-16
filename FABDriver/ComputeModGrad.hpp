//This class computes the modulus of the magnitude of chi
#ifndef COMPUTEMODGRAD_HPP_
#define COMPUTEMODGRAD_HPP_

#include "UserVariables.hpp"
#include "tensor.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Cell.hpp"

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
   void compute(Cell current_cell)
   {
       tensor<1,data_t> d1_arr[c_NUM];
       FOR1(idir) m_deriv.diff1(d1_arr, current_cell, idir);

       std::array<data_t, c_NUM> mod_d1_arr = {0.};
       FORVARS(ivar)
       {
           FOR1(idir)
           {
               mod_d1_arr[ivar] += d1_arr[ivar][idir]*d1_arr[ivar][idir];
           }
           mod_d1_arr[ivar] = sqrt(mod_d1_arr[ivar]);
       }

       // Write back into the flattened Chombo box
       m_driver.store_vars(mod_d1_arr, current_cell);
   }

};

#endif /* COMPUTEMODGRAD_HPP_ */
