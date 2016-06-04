#include "FABDriver.hpp"
#include "simd.hpp"

//template <typename... param_types>
//FABDriver::FABDriver(param_types... params) : m_compute(compute_t(std::forward<param_types>(params)..., *this)) {};

template <class compute_t>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out)
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
            m_compute.template compute<simd<double> >(x, y, z);
         }

         // REMAINDER LOOP
#pragma novector
         for (int x = x_simd_max + simd<double>::simd_len; x <= m_out_hi[0]; ++x)
         {
            m_compute.template compute<double>(x, y, z);
         }
      }
}

