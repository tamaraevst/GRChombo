#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

//This class offers a nan-check including debugging information when it happens

#include "user_enum.hpp"
#include "FABDriverBase.hpp"

class NanCheck
{
protected:
   const FABDriverBase& m_driver;

public:
   NanCheck(const FABDriverBase& driver) :
      m_driver (driver)
   {}

   void compute(int ix, int iy, int iz)
   {
      idx_t<double> idx = m_driver.in_idx(ix, iy, iz);
      FORVARS(i)
      {
          const double val = m_driver.local_vars(idx,i);
          if ( isnan(val) )
              MayDay::Error("values have become nan");
      }
   }
};

#endif /* NANCHECK_HPP_ */
