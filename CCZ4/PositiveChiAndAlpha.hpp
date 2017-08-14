//This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "UserVariables.hpp"
#include "MiscUtils.hpp"
#include "simd.hpp"
#include "Cell.hpp"

class PositiveChiAndAlpha
{
public:
   template <class data_t>
   void compute(Cell<data_t> current_cell)
   {
      auto chi = current_cell.load_vars(c_chi);
      auto lapse = current_cell.load_vars(c_lapse);

      MIN_CUT_OFF(chi, 1e-4);
      MIN_CUT_OFF(lapse, 1e-4);

      current_cell.store_vars(chi, c_chi);
      current_cell.store_vars(lapse, c_lapse);
   }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
