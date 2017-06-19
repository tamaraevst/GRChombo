#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

//This class offers a nan-check including debugging information when it happens

#include "UserVariables.hpp"
#include "Cell.hpp"

class NanCheck
{
protected:
   const std::string m_error_info = "NanCheck";
   const double m_max_abs = 1e20;

public:
   NanCheck() {}

   NanCheck(const std::string a_error_info) : m_error_info (a_error_info) {}

   NanCheck(const std::string a_error_info, const double a_max_abs) :
       m_error_info (a_error_info), m_max_abs (a_max_abs) {}

   void compute( Cell<double> current_cell )
   {
      bool stop = false;
      FORVARS(i)
      {
          double val;
          current_cell.local_vars(val,i);
          if ( std::isnan(val) || abs(val) > m_max_abs) stop = true;
      }
      if (stop)
      {
#pragma omp single
          {
              pout() << m_error_info << "::Values have become nan. The current state is: " << endl;
              FORVARS(i)
              {
                  pout() << UserVariables::variable_names[i] << ": " << current_cell.local_vars(i) << endl;
              }
              pout() << "Integer coordinates: " << current_cell.get_int_vect() << endl;
          }
          MayDay::Error("Values have become nan.");
      }
   }
};

#endif /* NANCHECK_HPP_ */
