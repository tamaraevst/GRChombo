#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

//This class offers a nan-check including debugging information when it happens

#include "UserVariables.hpp"
#include "FABDriverBase.hpp"

class NanCheck
{
protected:
   const FABDriverBase& m_driver;
   const std::string m_error_info = "NanCheck";
   const double m_max_abs = 1e20;

public:
   NanCheck(const FABDriverBase& driver) : m_driver (driver) {}

   NanCheck(const FABDriverBase& driver, const std::string a_error_info) : m_driver (driver), m_error_info (a_error_info) {}

   NanCheck(const FABDriverBase& driver, const std::string a_error_info, const double a_max_abs) :
       m_driver (driver), m_error_info (a_error_info), m_max_abs (a_max_abs) {}

   void compute(int ix, int iy, int iz)
   {
      idx_t<double> idx = m_driver.in_idx(ix, iy, iz);
      bool stop = false;
      FORVARS(i)
      {
          const double val = m_driver.local_vars(idx,i);
          if ( std::isnan(val) || abs(val) > m_max_abs) stop = true;
      }
      if (stop)
      {
#pragma omp single
          {
              pout() << m_error_info << "::Values have become nan. The current state is: " << endl;
              FORVARS(i)
              {
                  pout() << UserVariables::variable_names[i] << ": " << m_driver.local_vars(idx,i) << endl;
              }
              pout() << "ix: " << ix << " iy: " << iy << " iz: " << iz << endl;
          }
          MayDay::Error("Values have become nan.");
      }
   }
};

#endif /* NANCHECK_HPP_ */
