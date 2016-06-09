//This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "user_enum.hpp"
#include "simd.hpp"
#include "FABDriverBase.hpp"

class PositiveChiAndAlpha
{
protected:
   const FABDriverBase& m_driver;

public:
   PositiveChiAndAlpha(const FABDriverBase& driver) :
      m_driver (driver)
   {}

   template <class data_t>
   void compute(int x, int y, int z)
   {
      idx_t<data_t> idx = m_driver.in_idx(x, y, z);

      auto vars_arr = m_driver.local_vars(idx);

      auto chi_is_too_small = simd_compare_lt(vars_arr[c_chi], 1e-4);
      vars_arr[c_chi] = simd_conditional(chi_is_too_small, 1e-4, vars_arr[c_chi]);

      auto lapse_is_too_small = simd_compare_lt(vars_arr[c_lapse], 1e-4);
      vars_arr[c_lapse] = simd_conditional(lapse_is_too_small, 1e-4, vars_arr[c_lapse]);

      idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_chi])[out_idx]   = vars_arr[c_chi];
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_lapse])[out_idx] = vars_arr[c_lapse];
   }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
