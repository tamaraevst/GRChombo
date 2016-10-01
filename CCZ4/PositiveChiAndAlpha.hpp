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
      data_t chi, lapse;
      m_driver.local_vars(chi, c_chi);
      m_driver.local_vars(lapse, c_lapse);

      auto chi_is_too_small = simd_compare_lt(chi, 1e-4);
      chi = simd_conditional(chi_is_too_small, 1e-4, chi);

      auto lapse_is_too_small = simd_compare_lt(lapse, 1e-4);
      lapse = simd_conditional(lapse_is_too_small, 1e-4, lapse);

      m_driver.store_vars(chi, c_chi);
      m_driver.store_vars(lapse, c_lapse);
   }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
