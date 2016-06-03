#include "FABDriverBase.hpp"
#include "simd.hpp"

template <class data_t>
void
FABDriverBase::local_vars(int idx, data_t (&out)[c_NUM]) const const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        out[i] = SIMDIFY<data_t>(m_in_ptr[i])[idx];
    }
}

template <class data_t>
void
FABDriverBase::diff1(int idx, int stride, double dx, data_t (&out)[c_NUM]) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far  = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        out[i] = (  weight_far  * in[idx - 2*stride]
                     - weight_near * in[idx - stride]
                     + weight_near * in[idx + stride]
                     - weight_far  * in[idx + 2*stride]) / dx;
    }
}

template <class data_t>
void
FABDriverBase::diff2(int idx, int stride, double dx, data_t (&out)[c_NUM]) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far   = 8.33333333333333333333e-2;
        data_t weight_near  = 1.33333333333333333333e+0;
        data_t weight_local = 2.50000000000000000000e+0;

        out[i] = (- weight_far   * in[idx - 2*stride]
                     + weight_near  * in[idx - stride]
                     - weight_local * in[idx]
                     + weight_near  * in[idx + stride]
                     - weight_far   * in[idx + 2*stride]) / (dx*dx);
    }
}

template <class data_t>
void
FABDriverBase::mixed_diff2(int idx, int stride1, int stride2, double dx, data_t (&out)[c_NUM]) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

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
                     + weight_far_far   * in[idx + 2*stride1 + 2*stride2]) / (dx*dx);
    }
}

template <class data_t>
void
FABDriverBase::advection(int idx, double dx, const tensor<1, data_t>& shift, data_t (&out)[c_NUM]) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        out[i] = 0;
    }

    for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
    {
        const int stride = m_stride[dim];
        const auto shift_val = shift[dim];
        const auto shift_positive = simd_compare_gt(shift_val, 0.0);

        for (int i = 0; i < c_NUM; ++i)
        {
            const auto in = SIMDIFY<data_t>(m_in_ptr[i]);
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
                                  + weight_4 * in[idx + 3*stride]) / dx;

            data_t downwind;
            downwind = shift_val * (- weight_4 * in[idx - 3*stride]
                                    - weight_3 * in[idx - 2*stride]
                                    - weight_2 * in_left
                                    - weight_1 * in_centre
                                    - weight_0 * in_right) / dx;

            out[i] += simd_conditional(shift_positive, upwind , downwind);
        }
    }
}

template <class data_t>
void
FABDriverBase::dissipation(int idx, double dx, data_t (&out)[c_NUM]) const
{
    for (int i = 0; i < c_NUM; ++i)
    {
        out[i] = 0;
    }

    for (int dim = 0; dim < IDX_SPACEDIM; ++dim)
    {
        const int stride = m_stride[dim];
        for (int i = 0; i < c_NUM; ++i)
        {
            auto in = SIMDIFY<data_t>(m_in_ptr[i]);

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
                          + weight_vfar  * in[idx + 3*stride]) / dx;
        }
    }
}
