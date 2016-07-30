#ifndef VARSBASE_HPP_
#define VARSBASE_HPP_

//This class is a base class for variables type in compute classes.
//Inherit from it to get the functionality of loading grid variables into the variable type
//Derivatives are implemented using tensor<N, data_t>

#include "StackVector.hpp"
#include "IndexApplicator.hpp"

template <class var_t>
class VarsBase
{
protected:
     //holds pointers that can be used for assignment.
     //Every component has 2 pointers (in case we have symmetric tensors)
     StackVector<var_t*,2> m_assignment_ptrs[c_NUM];

public:
     //Writes data directly into the variable corresponding to ivar
     //if this variables has multiple components (e.g. if it is an array of derivatives)
     //the data can be written directly into these components by specifying
     //an arbitrary number of icomps
     template <class data_t, typename... Ints>
     ALWAYS_INLINE
     void assign(const data_t& data, int ivar, Ints&&... icomps)
     {
         for (int i=0; i<m_assignment_ptrs[ivar].get_ncomp(); ++i)
         {
             IndexApplicator<Ints...>::apply(
                 *m_assignment_ptrs[ivar][i], std::forward<Ints>(icomps)...) = data;
         }
     }

     //Adds the input data, otherwise the behaviour is similar to assign
     template <class data_t, typename... Ints>
     ALWAYS_INLINE
     void plus(const data_t& data, int ivar, Ints&&... icomps)
     {
         for (int i=0; i<m_assignment_ptrs[ivar].get_ncomp(); ++i)
         {
             IndexApplicator<Ints...>::apply(
                 *m_assignment_ptrs[ivar][i], std::forward<Ints>(icomps)...) += data;
         }
     }

     //Assigns data_t to each component of each variable
     template <class data_t>
     ALWAYS_INLINE
     void assign(const data_t& data)
     {
         FORVARS(ivar) assign(data, ivar);
     }

     //This functions copies data from an input array
     //Use it only if you must - wherever possible write the data in directly
     template <class data_t>
     ALWAYS_INLINE
     void
     assign(const var_t (&data)[c_NUM])
     {
         FORVARS(ivar) assign(data[ivar], ivar);
     }


};
#endif /* VARSBASE_HPP_ */
