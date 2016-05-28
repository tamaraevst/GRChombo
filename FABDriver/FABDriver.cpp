
template <typename... param_types>
FABDriver::FABDriver(param_types... params) :
    m_compute(std::forward<param_types>(params)...)
{}

    void execute(const FArrayBox& in, FArrayBox& out)
    {
        // dataPtr in Chombo does CH_assert bound check
        // which we don't want to do in a loop
        for (int i = 0; i < c_NUM; ++i)
        {
            m_in_ptr[i] = in.dataPtr(i);
            m_out_ptr[i] = out.dataPtr(i);
        }

        m_in_lo = in.loVect();
        m_in_hi = in.hiVect();
        m_stride[0] = 1;
        m_stride[1] = m_in_hi[0]-m_in_lo[0]+1;
    #if CH_SPACEDIM >= 3
        m_stride[2] = (m_in_hi[1]-m_in_lo[1]+1)*m_stride[1];
    #endif
    #if CH_SPACEDIM >= 4
        #error "TODO: Implement CH_SPACEDIM >= 4"
    #endif

        m_out_lo = out.loVect();
        m_out_hi = out.hiVect(); 
        m_out_stride[0] = 1;
        m_out_stride[1] = m_out_hi[0]-m_out_lo[0]+1;
    #if CH_SPACEDIM >= 3
        m_out_stride[2] = (m_out_hi[1]-m_out_lo[1]+1)*m_out_stride[1];
    #endif
    #if CH_SPACEDIM >= 4
        #error "TODO: Implement CH_SPACEDIM >= 4"
    #endif

    #pragma omp parallel for default(shared) collapse(CH_SPACEDIM-1)
    #if CH_SPACEDIM >= 4
        #error "TODO: Implement CH_SPACEDIM >= 4"
    #endif
    #if CH_SPACEDIM >= 3
        for (int z = m_out_lo[2]; z <= m_out_hi[2]; ++z)
    #endif
        for (int y = m_out_lo[1]; y <= m_out_hi[1]; ++y)
        {
            int x_simd_max = m_out_lo[0] + simd<double>::simd_len * (((m_out_hi[0] - m_out_lo[0] + 1) / simd<double>::simd_len) - 1);

            // SIMD LOOP
            #pragma novector
            for (int x = m_out_lo[0]; x <= x_simd_max; x += simd<double>::simd_len)
            {
                m_compute.compute<simd<double> >(x, y, z);
            }

            // REMAINDER LOOP
            #pragma novector
            for (int x = x_simd_max + simd<double>::simd_len; x <= m_out_hi[0]; ++x)
            {
                m_compute.compute<double>(x, y, z);
            }
        }
    }

    template <class data_t>
void
CCZ4::local_vars(int idx, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = SIMDIFY<data_t>(m_in_ptr[i])[idx];
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::diff1(int idx, int stride, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far  = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        result[i] = (  weight_far  * in[idx - 2*stride]
                     - weight_near * in[idx - stride]
                     + weight_near * in[idx + stride]
                     - weight_far  * in[idx + 2*stride]) / m_dx;
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::diff2(int idx, int stride, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far   = 8.33333333333333333333e-2;
        data_t weight_near  = 1.33333333333333333333e+0;
        data_t weight_local = 2.50000000000000000000e+0;

        result[i] = (- weight_far   * in[idx - 2*stride]
                     + weight_near  * in[idx - stride]
                     - weight_local * in[idx]
                     + weight_near  * in[idx + stride]
                     - weight_far   * in[idx + 2*stride]) / (m_dx*m_dx);
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::mixed_diff2(int idx, int stride1, int stride2, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        auto in = SIMDIFY<data_t>(m_in_ptr[i]);

        data_t weight_far_far   = 6.94444444444444444444e-3;
        data_t weight_near_far  = 5.55555555555555555556e-2;
        data_t weight_near_near = 4.44444444444444444444e-1;

        result[i] = (  weight_far_far   * in[idx - 2*stride1 - 2*stride2]
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
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::advection(int idx, const tensor<1, data_t>& shift, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = 0;
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
                                  + weight_4 * in[idx + 3*stride]) / m_dx;

            data_t downwind;
            downwind = shift_val * (- weight_4 * in[idx - 3*stride]
                                    - weight_3 * in[idx - 2*stride]
                                    - weight_2 * in_left
                                    - weight_1 * in_centre
                                    - weight_0 * in_right) / m_dx;

            result[i] += simd_conditional(shift_positive, upwind , downwind);
        }
    }
    demarshall(result, out);
}

template <class data_t>
void
CCZ4::dissipation(int idx, vars_t<data_t>& out)
{
    data_t result[c_NUM];
    for (int i = 0; i < c_NUM; ++i)
    {
        result[i] = 0;
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

            result[i] += ( weight_vfar   * in[idx - 3*stride]
                          - weight_far   * in[idx - 2*stride]
                          + weight_near  * in[idx - stride]
                          - weight_local * in[idx]
                          + weight_near  * in[idx + stride]
                          - weight_far   * in[idx + 2*stride]
                          + weight_vfar  * in[idx + 3*stride]) / m_dx;
        }
    }
    demarshall(result, out);
}
}