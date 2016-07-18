#ifndef VARSBASE_HPP_
#define VARSBASE_HPP_

//This class is a base class for variables type in compute classes.
//Inherit from it to get the functionality of loading grid variables into the variable type
//Derivatives are implemented using tensor<N, data_t>

#include "StackVector.hpp"

template <class data_t>
class VarsBase
{
protected:
     //holds pointers that can be used for assignement.
     //Every component has 2 pointers (in case we have symmetric tensors)
     StackVector<data_t*,2> m_assignement_ptrs[c_NUM];

     void
     assign(data_t& data, int icomp)
     {
         for (int i=0; i<m_assignement_ptrs[icomp].get_ncomp(); ++i)
         {
             m_assignement_ptrs[icomp][i] = data;
         }
     }

     void
     assign(data_t& data[c_NUM])
     {
         FORCOMPS(icomp) assign(data[icomp], icomp);
     }
};
#endif /* VARSBASE_HPP_ */
