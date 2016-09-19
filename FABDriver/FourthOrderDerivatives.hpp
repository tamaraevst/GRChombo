#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include "VarsBase.hpp"
#include "FABDriverBase.hpp"
#include <array>

class FourthOrderDerivatives
{
public:
    const double m_dx;
protected:
    const FABDriverBase& m_driver;

public:
    FourthOrderDerivatives(double dx, const FABDriverBase& driver) :
        m_dx (dx),
        m_driver (driver)
    {}

    template <class data_t>
    ALWAYS_INLINE
    data_t
    diff1(const double *in_ptr, idx_t<data_t> idx, const int stride) const
    {
        auto in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far  = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        //NOTE: if you have been sent here by the debugger because of EXC_BAD_ACCESS
        //or something similar you might be trying to take derivatives without ghost points.
        return (  weight_far  * in[idx - 2*stride]
                    - weight_near * in[idx - stride]
                    + weight_near * in[idx + stride]
                    - weight_far  * in[idx + 2*stride]) / m_dx;
    }

    //Writes directly into the vars object - use this wherever possible
    template <class data_t>
    void
    diff1(VarsBase< tensor<1, data_t> >& d1, idx_t<data_t> idx, int direction) const
    {
        int stride = m_driver.m_in_stride[direction];
        FORVARS(i)
        {
            d1.assign(diff1(m_driver.m_in_ptr[i], idx, stride), i, direction);
        }
    }

    template <class data_t>
    void
    diff1(tensor<1,data_t> (&diffArray)[c_NUM], idx_t<data_t> idx, int direction) const
    {
        int stride = m_driver.m_in_stride[direction];
        FORVARS(i)
        {
            diffArray[i][direction] = diff1(m_driver.m_in_ptr[i],idx,stride);
        }
    }

    template <class data_t>
    ALWAYS_INLINE
    data_t
    diff2(const double *in_ptr, idx_t<data_t> idx, const int stride) const
    {
            auto in = SIMDIFY<data_t>(in_ptr);

            data_t weight_far   = 8.33333333333333333333e-2;
            data_t weight_near  = 1.33333333333333333333e+0;
            data_t weight_local = 2.50000000000000000000e+0;

            return   (- weight_far   * in[idx - 2*stride]
                         + weight_near  * in[idx - stride]
                         - weight_local * in[idx]
                         + weight_near  * in[idx + stride]
                         - weight_far   * in[idx + 2*stride]) / (m_dx*m_dx);
    }

    //Writes 2nd deriv directly into the vars object - use this wherever possible
    template <class data_t>
    void
    diff2(VarsBase< tensor<2,data_t> >& d2, idx_t<data_t> idx, int direction) const
    {
        int stride = m_driver.m_in_stride[direction];
        FORVARS(i)
        {
            d2.assign(diff2(m_driver.m_in_ptr[i], idx, stride), i, direction, direction);
        }
    }

    template <class data_t>
    void
    diff2(tensor<2,data_t> (&diffArray)[c_NUM], idx_t<data_t> idx, int direction) const
    {
        int stride = m_driver.m_in_stride[direction];
        FORVARS(i)
        {
            diffArray[i][direction][direction] = diff2(m_driver.m_in_ptr[i],idx,stride);
        }
    }

    template <class data_t>
    ALWAYS_INLINE
    data_t
    mixed_diff2(const double *in_ptr, idx_t<data_t> idx, const int stride1, const int stride2) const
    {
            auto in = SIMDIFY<data_t>(in_ptr);

            data_t weight_far_far   = 6.94444444444444444444e-3;
            data_t weight_near_far  = 5.55555555555555555556e-2;
            data_t weight_near_near = 4.44444444444444444444e-1;

            return (  weight_far_far   * in[idx - 2*stride1 - 2*stride2]
                         - weight_near_far  * in[idx - 2*stride1 - stride2]
                         + weight_near_far  * in[idx - 2*stride1 + stride2]
                         - weight_far_far   * in[idx - 2*stride1 + 2*stride2]

                         - weight_near_far  * in[idx - stride1 - 2*stride2]
                         + weight_near_near * in[idx - stride1 - stride2]
                         - weight_near_near * in[idx - stride1 + stride2]
                         + weight_near_far  * in[idx - stride1 + 2*stride2]

                         + weight_near_far  * in[idx + stride1 - 2*stride2]
                         - weight_near_near * in[idx + stride1 - stride2]
                         + weight_near_near * in[idx + stride1 + stride2]
                         - weight_near_far  * in[idx + stride1 + 2*stride2]

                         - weight_far_far   * in[idx + 2*stride1 - 2*stride2]
                         + weight_near_far  * in[idx + 2*stride1 - stride2]
                         - weight_near_far  * in[idx + 2*stride1 + stride2]
                         + weight_far_far   * in[idx + 2*stride1 + 2*stride2]) / (m_dx*m_dx);
    }

    template <class data_t>
    void
    mixed_diff2(VarsBase< tensor<2,data_t> >& d2, idx_t<data_t> idx, int direction1, int direction2) const
    {
        int stride1 = m_driver.m_in_stride[direction1];
        int stride2 = m_driver.m_in_stride[direction2];
        FORVARS(i)
        {
            data_t tmp = mixed_diff2(m_driver.m_in_ptr[i], idx, stride1, stride2);
            d2.assign(tmp, i, direction1, direction2);
            d2.assign(tmp, i, direction2, direction1);
        }
    }

    template <class data_t>
    void
    mixed_diff2(tensor<2,data_t> (&diffArray)[c_NUM], idx_t<data_t> idx, int direction1, int direction2) const
    {
        int stride1 = m_driver.m_in_stride[direction1];
        int stride2 = m_driver.m_in_stride[direction2];
        FORVARS(i)
        {
            data_t diff2_value = mixed_diff2(m_driver.m_in_ptr[i],idx,stride1, stride2);
            diffArray[i][direction1][direction2] = diff2_value;
            diffArray[i][direction2][direction1] = diff2_value;
        }
    }

protected: //Let's keep this protected ... we may want to change the advection calculation
    template <class data_t>
    ALWAYS_INLINE
    data_t
    advection_term(const double *in_ptr, idx_t<data_t> idx, const data_t& vec_comp, const int stride, const data_t shift_positive) const
    {
        const auto in = SIMDIFY<data_t>(in_ptr);
        const data_t in_left = in[idx - stride];
        const data_t in_centre = in[idx];
        const data_t in_right = in[idx + stride];

        data_t weight_0 = -2.50000000000000000000e-1;
        data_t weight_1 = -8.33333333333333333333e-1;
        data_t weight_2 = +1.50000000000000000000e+0;
        data_t weight_3 = -5.00000000000000000000e-1;
        data_t weight_4 = +8.33333333333333333333e-2;

        data_t upwind;
        upwind = vec_comp * (  weight_0 * in_left
                          + weight_1 * in_centre
                          + weight_2 * in_right
                          + weight_3 * in[idx + 2*stride]
                          + weight_4 * in[idx + 3*stride]) / m_dx;

        data_t downwind;
        downwind = vec_comp * (- weight_4 * in[idx - 3*stride]
                          - weight_3 * in[idx - 2*stride]
                          - weight_2 * in_left
                          - weight_1 * in_centre
                          - weight_0 * in_right) / m_dx;

        return simd_conditional(shift_positive, upwind , downwind);
    }

public:
    template <class data_t>
    void
    add_advection(VarsBase<data_t>& vars, idx_t<data_t> idx, const data_t& vec_comp, const int dir) const
    {
        const int stride = m_driver.m_in_stride[dir];
        const data_t shift_positive = simd_compare_gt(vec_comp, 0.0);
        FORVARS(ivar)
        {
            vars.plus(advection_term(m_driver.m_in_ptr[ivar],idx, vec_comp, stride, shift_positive),ivar);
        }
    }

    template <class data_t>
    void
    add_advection(data_t (&out)[c_NUM], idx_t<data_t> idx, const data_t& vec_comp, const int dir) const
    {
        const int stride = m_driver.m_in_stride[dir];
        const data_t shift_positive = simd_compare_gt(vec_comp, 0.0);
        FORVARS(ivar)
        {
            out[ivar] += advection_term(m_driver.m_in_ptr[ivar],idx, vec_comp, stride, shift_positive);
        }
    }

    template <class data_t>
    ALWAYS_INLINE
    data_t
    dissipation_term(const double *in_ptr, idx_t<data_t> idx, const int stride) const
    {
        const auto in = SIMDIFY<data_t>(in_ptr);
        data_t weight_vfar  = 1.56250e-2;
        data_t weight_far   = 9.37500e-2;
        data_t weight_near  = 2.34375e-1;
        data_t weight_local = 3.12500e-1;

        return ( weight_vfar   * in[idx - 3*stride]
                 - weight_far   * in[idx - 2*stride]
                 + weight_near  * in[idx - stride]
                 - weight_local * in[idx]
                 + weight_near  * in[idx + stride]
                 - weight_far   * in[idx + 2*stride]
                 + weight_vfar  * in[idx + 3*stride]) / m_dx;
    }

    template <class data_t>
    void
    add_dissipation(VarsBase<data_t>& vars, const double factor, idx_t<data_t> idx, const int direction) const
    {
        const int stride = m_driver.m_in_stride[direction];
        FORVARS(ivar)
        {
            vars.plus(factor*dissipation_term(m_driver.m_in_ptr[ivar], idx, stride), ivar);
        }
    }

    template <class data_t>
    void
    add_dissipation(data_t (&out)[c_NUM], const double factor, idx_t<data_t> idx, const int direction) const
    {
        const int stride = m_driver.m_in_stride[direction];
        FORVARS(ivar)
        {
            out[ivar] += factor*dissipation_term(m_driver.m_in_ptr[ivar], idx, stride);
        }
    }
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
