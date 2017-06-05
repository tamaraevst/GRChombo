#ifndef PHIANDKTAGGINGCRITERION_HPP_
#define PHIANDKTAGGINGCRITERION_HPP_

class PhiAndKTaggingCriterion
{
protected:
   const double m_dx;
   const double m_threshold_phi;
   const double m_threshold_K;
   const FourthOrderDerivatives m_deriv;

public:
   PhiAndKTaggingCriterion(double dx, double threshold_phi, double threshold_K)
    : m_dx (dx), m_deriv (dx), m_threshold_phi (threshold_phi), m_threshold_K (threshold_K) {};

   template <class data_t>
   void compute(Cell current_cell)
   {
       auto phi = current_cell.template local_vars<data_t>(c_phi);
       tensor<1,data_t> d1_phi;
       FOR1(idir) m_deriv.diff1(d1_phi, current_cell, idir, c_phi);

       auto K = current_cell.template local_vars<data_t>(c_K);
       tensor<1,data_t> d1_K;
       FOR1(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

       data_t mod_d1_phi = 0;
       data_t mod_d1_K = 0;
       FOR1(idir)
       {
           mod_d1_phi += d1_phi[idir]*d1_phi[idir];
           mod_d1_K += d1_K[idir]*d1_K[idir];
       }

       data_t criterion = m_dx * (1.0/m_threshold_phi*(sqrt(mod_d1_phi)/pow(phi,2)) 
                              + 1.0/m_threshold_K*(sqrt(mod_d1_K)/pow(K,2)));

       // Write back into the flattened Chombo box
       current_cell.store_vars(criterion, 0);
   }
};

#endif /* PHIANDKTAGGINGCRITERION_HPP_ */
