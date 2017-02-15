#ifndef VARSBASE_HPP_
#define VARSBASE_HPP_

#include "StackVector.hpp"
#include "IndexApplicator.hpp"
#include "GRUtils.hpp"
#include "CH_assert.H"
#include "Interval.H"
#include "UserVariables.hpp"
#include "tensor.hpp"

///This class is a base class for variables types in compute classes.
/**Inherit from it to get the functionality of loading grid variables into your new variables class.
 * VarsBase contains a list of c_NUM (as defined in UserVariables.hpp) assignment pointers which can map a number in the enum
 * defined in UserVariables.hpp to a variable of the variables class that inherits from VarsBase.
 * These assignment pointers are then used to load values from the Chombo grid into the variables class.
 *
 * e.g. A child Vars of VarsBase may have a member variable 'chi'. VarsBase m_assignement_ptrs[c_chi] should then be
 * set up to point to this variable chi. To do this use the functions define_enum_mapping.
 *
 * NB: A VarsBase object cannot be copied as this almost certainly causes troubles with the assignement pointers
 * pointing to the old data. If you must copy a child of VarsBase you have to write an explicit copy constructor
 * and = operator for the child. These should create a new VarsBase object and explicitly fix the assignment pointers.
 */
template <class var_t>
class VarsBase
{
protected:
     //holds pointers that can be used for assignment.
     //Every component has 2 pointers (in case we have symmetric tensors)
     StackVector<var_t*,2> m_assignment_ptrs[c_NUM];

public:
     VarsBase(VarsBase const &) = delete; //!< VarsBase cannot be copied (see class description).
     void operator=(VarsBase const &vars) = delete; //!< VarsBase cannot be copied (see class description).
     VarsBase() {} //!< Default constructor. Leaves assignment pointers empty.

     //This function assigns a mapping from a chombo grid variable to a local
     //variable (i.e. a variable that only exists for the current cell).
     //This can be used up to two times (e.g. for
     //symmetric tensor h12 would map to both vars.h[0][1] and vars.h[1][0])
     void define_enum_mapping(int a_enum_component, var_t& a_var)
     {
         m_assignment_ptrs[a_enum_component].push_back(&a_var);
     }

     void define_enum_mapping(Interval a_enum_vector_components, tensor<1, var_t>& a_var)
     {
         CH_assert (a_enum_vector_components.size() == IDX_SPACEDIM);
         int icomp = 0;
         for (int i=a_enum_vector_components.begin(); i<=a_enum_vector_components.end(); ++i)
         {
             m_assignment_ptrs[i].push_back(&(a_var[icomp]));
             ++icomp;
         }
     }

     void define_symmetric_enum_mapping(Interval a_enum_tensor_components, tensor<2, var_t>& a_var)
     {
         CH_assert (a_enum_tensor_components.size() == IDX_SPACEDIM*(IDX_SPACEDIM+1)/2.);
         int start_comp = a_enum_tensor_components.begin();
#if IDX_SPACEDIM == 3
         m_assignment_ptrs[start_comp  ].push_back(&(a_var[0][0]));

         m_assignment_ptrs[start_comp+1].push_back(&(a_var[0][1]));
         m_assignment_ptrs[start_comp+1].push_back(&(a_var[1][0]));

         m_assignment_ptrs[start_comp+2].push_back(&(a_var[0][2]));
         m_assignment_ptrs[start_comp+2].push_back(&(a_var[2][0]));

         m_assignment_ptrs[start_comp+3].push_back(&(a_var[1][1]));

         m_assignment_ptrs[start_comp+4].push_back(&(a_var[1][2]));
         m_assignment_ptrs[start_comp+4].push_back(&(a_var[2][1]));

         m_assignment_ptrs[start_comp+5].push_back(&(a_var[2][2]));
#else
#error IDX_SPACEDIM not equal to three not implemented yet...
#endif
     }

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

     ALWAYS_INLINE
     bool variable_defined(int ivar) const
     {
         return m_assignment_ptrs[ivar].get_ncomp();
     }

     template <class data_t>
     ALWAYS_INLINE
     const data_t& read(int ivar) const
     {
         CH_assert(variable_defined(ivar));
         //the assignement pointer points to several variables with the same value e.g. symmetric matrix components.
         //The 0 component must exists so just take this one.
         return *m_assignment_ptrs[ivar][0];
     }
};
#endif /* VARSBASE_HPP_ */
