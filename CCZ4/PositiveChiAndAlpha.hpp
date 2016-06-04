//This compute class enforces the positive chi and alpha condition
#include "user_enum.hpp"
#include "simd.hpp"

class PositiveChiAndAlpha
{
   public:
      PositiveChiAndAlpha(const FABDriverBase& driver): m_driver (driver){};

      template <class data_t>
         void compute(int x, int y, int z)
         {
            const int idx = m_driver.m_stride[2]*(z-m_driver.m_in_lo[2]) + m_driver.m_stride[1]*(y-m_driver.m_in_lo[1]) + (x-m_driver.m_in_lo[0]);
            data_t varsArr[c_NUM];
            m_driver.local_vars(idx, varsArr);

            auto chiTooSmall = simd_compare_lt(varsArr[c_chi], 1e-4);
            varsArr[c_chi] = simd_conditional(chiTooSmall, 1e-4, varsArr[c_chi]);

            auto lapseTooSmall = simd_compare_lt(varsArr[c_lapse], 1e-4);
            varsArr[c_lapse] = simd_conditional(lapseTooSmall, 1e-4, varsArr[c_lapse]);
         }

   protected:
      const FABDriverBase& m_driver;
};
