#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "user_enum.hpp"
#include "tensor.hpp"
#include <array>

struct FourthOrderDerivatives
{
    const double m_dx;
    const FABDriverBase& m_driver;

    FourthOrderDerivatives(double dx, const FABDriverBase& driver) :
        m_dx (dx),
        m_driver (driver)
    {}

    template <class data_t>
    std::array<data_t, c_NUM>
    diff1(idx_t<data_t> idx, int direction) const
    {
        std::array<data_t, c_NUM> out;
        int stride = m_driver.m_in_stride[direction];

        for (int i = 0; i < c_NUM; ++i)
        {
            auto in = SIMDIFY<data_t>(m_driver.m_in_ptr[i]);

            data_t weight_far  = 8.33333333333333333333e-2;
            data_t weight_near = 6.66666666666666666667e-1;

            //NOTE: if you have been sent here by the debugger because of EXC_BAD_ACCESS
            //or something similar you might be trying to take derivatives without ghost points.
            out[i] = (  weight_far  * in[idx - 2*stride]
                         - weight_near * in[idx - stride]
                         + weight_near * in[idx + stride]
                         - weight_far  * in[idx + 2*stride]) / m_dx;
        }
        return out;
    }

    template <class data_t>
    std::array<data_t, c_NUM>
    diff2(idx_t<data_t> idx, int direction) const
    {
        std::array<data_t, c_NUM> out;
        int stride = m_driver.m_in_stride[direction];

        for (int i = 0; i < c_NUM; ++i)
        {
            auto in = SIMDIFY<data_t>(m_driver.m_in_ptr[i]);

            data_t weight_far   = 8.33333333333333333333e-2;
            data_t weight_near  = 1.33333333333333333333e+0;
            data_t weight_local = 2.50000000000000000000e+0;

            out[i] = (- weight_far   * in[idx - 2*stride]
                         + weight_near  * in[idx - stride]
                         - weight_local * in[idx]
                         + weight_near  * in[idx + stride]
                         - weight_far   * in[idx + 2*stride]) / (m_dx*m_dx);
        }
        return out;
    }

    template <class data_t>
    std::array<data_t, c_NUM>
    mixed_diff2(idx_t<data_t> idx, int direction1, int direction2) const
    {
        std::array<data_t, c_NUM> out;
        int stride1 = m_driver.m_in_stride[direction1];
        int stride2 = m_driver.m_in_stride[direction2];

        for (int i = 0; i < c_NUM; ++i)
        {
            auto in = SIMDIFY<data_t>(m_driver.m_in_ptr[i]);

            data_t weight_far_far   = 6.94444444444444444444e-3;
            data_t weight_near_far  = 5.55555555555555555556e-2;
            data_t weight_near_near = 4.44444444444444444444e-1;

            out[i] = (  weight_far_far   * in[idx - 2*stride1 - 2*stride2]
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
        return out;
    }

    template <class data_t>
    std::array<data_t, c_NUM>
    advection(idx_t<data_t> idx, tensor<1, data_t> vec) const
    {
        std::array<data_t, c_NUM> out = { };

        for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
        {
            const int stride = m_driver.m_in_stride[dim];
            const auto shift_val = vec[dim];
            const auto shift_positive = simd_compare_gt(shift_val, 0.0);

            for (int i = 0; i < c_NUM; ++i)
            {
                const auto in = SIMDIFY<data_t>(m_driver.m_in_ptr[i]);
                const data_t in_left = in[idx - stride];
                const data_t in_centre = in[idx];
                const data_t in_right = in[idx + stride];

                data_t weight_0 = -2.50000000000000000000e-1;
                data_t weight_1 = -8.33333333333333333333e-1;
                data_t weight_2 = +1.50000000000000000000e+0;
                data_t weight_3 = -5.00000000000000000000e-1;
                data_t weight_4 = +8.33333333333333333333e-2;

                data_t upwind;
                upwind = shift_val * (  weight_0 * in_left
                                      + weight_1 * in_centre
                                      + weight_2 * in_right
                                      + weight_3 * in[idx + 2*stride]
                                      + weight_4 * in[idx + 3*stride]) / m_dx;

                data_t downwind;
                downwind = shift_val * (- weight_4 * in[idx - 3*stride]
                                        - weight_3 * in[idx - 2*stride]
                                        - weight_2 * in_left
                                        - weight_1 * in_centre
                                        - weight_0 * in_right) / m_dx;

                out[i] += simd_conditional(shift_positive, upwind , downwind);
            }
        }

        return out;
    }

    template <class data_t>
    std::array<data_t, c_NUM>
    dissipation(idx_t<data_t> idx) const
    {
        std::array<data_t, c_NUM> out = { };

        for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
        {
            const int stride = m_driver.m_in_stride[dim];
            for (int i = 0; i < c_NUM; ++i)
            {
                auto in = SIMDIFY<data_t>(m_driver.m_in_ptr[i]);

                data_t weight_vfar  = 1.56250e-2;
                data_t weight_far   = 9.37500e-2;
                data_t weight_near  = 2.34375e-1;
                data_t weight_local = 3.12500e-1;

                out[i] += ( weight_vfar   * in[idx - 3*stride]
                              - weight_far   * in[idx - 2*stride]
                              + weight_near  * in[idx - stride]
                              - weight_local * in[idx]
                              + weight_near  * in[idx + stride]
                              - weight_far   * in[idx + 2*stride]
                              + weight_vfar  * in[idx + 3*stride]) / m_dx;
            }
        }

        return out;
    }
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
