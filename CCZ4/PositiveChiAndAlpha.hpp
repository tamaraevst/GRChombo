//This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "UserVariables.hpp"
#include "simd.hpp"
#include "Cell.hpp"

class PositiveChiAndAlpha
{
public:
   template <class data_t>
   void compute(Cell current_cell)
   {
      data_t chi, lapse;
      current_cell.local_vars(chi, c_chi);
      current_cell.local_vars(lapse, c_lapse);

      auto chi_is_too_small = simd_compare_lt(chi, 1e-4);
      chi = simd_conditional(chi_is_too_small, 1e-4, chi);

      auto lapse_is_too_small = simd_compare_lt(lapse, 1e-4);
      lapse = simd_conditional(lapse_is_too_small, 1e-4, lapse);

      current_cell.store_vars(chi, c_chi);
      current_cell.store_vars(lapse, c_lapse);
   }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
