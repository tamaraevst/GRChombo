#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

//This class offers a nan-check including debugging information when it happens

#include "UserVariables.hpp"
#include "FABDriverBase.hpp"
#include "Cell.hpp"

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

   void compute( Cell current_cell )
   {
      bool stop = false;
      FORVARS(i)
      {
          double val;
          m_driver.local_vars(val,current_cell,i);
          if ( std::isnan(val) || abs(val) > m_max_abs) stop = true;
      }
      if (stop)
      {
#pragma omp single
          {
              pout() << m_error_info << "::Values have become nan. The current state is: " << endl;
              FORVARS(i)
              {
                  double val = m_driver.local_vars<double>(current_cell, i);
                  if ( std::isnan(val) || abs(val) > m_max_abs) pout() << "---> ";
                  else if ( abs(val) > m_max_abs/1e2 ) pout() << "   > ";
                  else pout() << "     ";
                  pout() << UserVariables::variable_names[i] << ": " << val << endl;
              }
              pout() << "Integer coordinates: " << current_cell.get_int_vect() << endl;
          }
          MayDay::Error("Values have become nan.");
      }
   }
};

#endif /* NANCHECK_HPP_ */
