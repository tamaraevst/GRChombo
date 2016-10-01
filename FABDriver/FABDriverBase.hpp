#ifndef FABDRIVERBASE_HPP_
#define FABDRIVERBASE_HPP_
//FABDriver and FABDriverBase implement looping of all points inside an FArrayBox
//with vectorisation and OpenMP.
//
//The reason why FABDriverBase is necessary is that if we want to inherit from one of the compute classes
//then the templating of the FABDriver with class compute_t becomes a problem when constructing the objects.
//FABDriverBase gives us all we need on the side of the computer but isn't templated.
//TODO: actually FABDriver doesn't need to be template over compute_t, only it's member functions

//The following file has to be provided by the user.
//It must include c_NUM, the total number of simulation components.
#include "UserVariables.hpp"
#include "always_inline.hpp"
#include "VarsBase.hpp"

#include <array>

class FABDriverBase
{
public: //TODO: these shouldn't be public ...
    const double *m_in_ptr[c_NUM];
    const int *m_in_lo;
    const int *m_in_hi;
    int m_in_stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];

protected:
    int m_in_idx;
    int m_out_idx;

public:
    //Some setters:
    ALWAYS_INLINE
    void set_idx(int ix, int iy, int iz);

    //Some getters:
    ALWAYS_INLINE
    int get_in_idx() const;

    ALWAYS_INLINE
    int get_out_idx() const;

    //The rest:
    template <class data_t>
    ALWAYS_INLINE
    data_t local_vars(int icomp) const;

    template <class data_t>
    ALWAYS_INLINE
    void local_vars(data_t& out, int icomp) const;

    template <class data_t>
    void local_vars(data_t (&out)[c_NUM]) const;

    template <class data_t>
    void local_vars(VarsBase<data_t>& vars) const;

    template <class data_t>
    void store_vars(const data_t& value, const int icomp) const;

    template <class data_t>
    void store_vars(const std::array<data_t, c_NUM>& values) const;

    template <class data_t>
    void store_vars(const VarsBase<data_t>& vars, const Interval a_comps) const;

    template <class data_t>
    void store_vars(const VarsBase<data_t>& vars) const;
};

#include "FABDriverBase.impl.hpp"

#endif /* FABDRIVERBASE_HPP_ */
