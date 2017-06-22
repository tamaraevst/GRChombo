#ifndef CHITAGGINGCRITERION_HPP_
#define CHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"

class ChiTaggingCriterion
{
protected:
   const double m_dx;
   const FourthOrderDerivatives m_deriv;
public:
   ChiTaggingCriterion(double dx) : m_dx (dx), m_deriv (dx) {};

   template <class data_t>
   void compute(Cell<data_t> current_cell)
   {
       auto chi = current_cell.local_vars(c_chi);
       tensor<1,data_t> d1_chi;
       FOR1(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

       data_t mod_d1_chi = 0;
       FOR1(idir)
       {
           mod_d1_chi += d1_chi[idir]*d1_chi[idir];
       }
       data_t criterion = m_dx * sqrt(mod_d1_chi)/pow(chi,2);

       // Write back into the flattened Chombo box
       current_cell.store_vars(criterion, 0);
   }
};

#endif /* CHITAGGINGCRITERION_HPP_ */