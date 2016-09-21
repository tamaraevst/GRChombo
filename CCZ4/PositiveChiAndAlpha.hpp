//This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "UserVariables.hpp"
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

      auto chi   = m_driver.local_vars(idx, c_chi);
      auto lapse = m_driver.local_vars(idx, c_lapse);

      auto chi_is_too_small = simd_compare_lt(chi, 1e-4);
      chi = simd_conditional(chi_is_too_small, 1e-4, chi);

      auto lapse_is_too_small = simd_compare_lt(lapse, 1e-4);
      lapse = simd_conditional(lapse_is_too_small, 1e-4, lapse);

      idx_t<data_t> out_idx = m_driver.out_idx(x, y, z);
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_chi])[out_idx]   = chi;
      SIMDIFY<data_t>(m_driver.m_out_ptr[c_lapse])[out_idx] = lapse;
   }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
