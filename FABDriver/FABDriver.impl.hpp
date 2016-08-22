#if !defined(FABDRIVER_HPP_)
#error "This file should only be included through FABDriver.hpp"
#endif

#ifndef FABDRIVER_IMPL_HPP_
#define FABDRIVER_IMPL_HPP_

#include "simd.hpp"

template <class compute_t>
template <typename... param_types>
FABDriver<compute_t>::FABDriver(param_types... params) : m_compute(*this, std::forward<param_types>(params)...) {}

//The following helpers allow the FABDriver to vectorise only if the provided compute class supports vectorisation.
//(MK: I wrote this so that the user doesn't have to worry about which compute classes can vectorise and which can't)
namespace SimdChecker
{
#warning MK: this is truly horrible code ... the calls are completely unreadable better to have two functions for the inner loop and switch them on or off
    template<class T, class data_t, typename... Args>
    static auto check_simd_compute(int, T& a_compute, Args... args) //passing an int as first argument makes the compiler try this version first
        -> decltype(std::declval<T>().template compute<data_t>(std::declval<Args>()...), std::true_type())
    {
        a_compute.template compute<simd<double> >(std::forward<Args>(args)...);
        return std::true_type();
    }
    template<class T, class, typename... Args>
    static auto check_simd_compute(char, T& a_compute, Args... args) //SFINAE fallback
        -> std::false_type
    {
        a_compute.compute(std::forward<Args>(args)...);
        return std::false_type();
    }

    template<class compute_t, class data_t>
    //returns std:true_type() if compute_t accepts simd<data_t> as template argument
    struct has_simd_compute :
        decltype(SimdChecker::check_simd_compute<compute_t, simd<data_t>, int, int, int> (0,std::declval<compute_t&>(),0,0,0)){};

    template<class data_t, class compute_t, typename... Args>
    void call_compute(compute_t& a_compute, Args... args)
    {
        check_simd_compute<compute_t,data_t>(0, a_compute, std::forward<Args>(args)... );
    }
}

template <class compute_t>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out, const Box& loop_box)
{
   //Makes sure we are not requesting data outside the box of 'out'
   CH_assert(out.box().contains(loop_box));

   // dataPtr in Chombo does CH_assert bound check
   // which we don't want to do in a loop
   for (int i = 0; i < c_NUM; ++i)
   {
      m_in_ptr[i] = in.dataPtr(i);
      m_out_ptr[i] = out.dataPtr(i);
   }

   m_in_lo = in.loVect();
   m_in_hi = in.hiVect();
   m_in_stride[0] = 1;
   m_in_stride[1] = m_in_hi[0]-m_in_lo[0]+1;
#if CH_SPACEDIM >= 3
   m_in_stride[2] = (m_in_hi[1]-m_in_lo[1]+1)*m_in_stride[1];
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

   const int* loop_lo = loop_box.loVect();
   const int* loop_hi = loop_box.hiVect();

#pragma omp parallel for default(shared) collapse(CH_SPACEDIM-1)
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif
#if CH_SPACEDIM >= 3
   for (int z = loop_lo[2]; z <= loop_hi[2]; ++z)
#endif
      for (int y = loop_lo[1]; y <= loop_hi[1]; ++y)
      {
         if (SimdChecker::has_simd_compute<compute_t, double>()) //If compute_t supports simd do a simd loop
         {
             int x_simd_max = loop_lo[0] + simd<double>::simd_len * (((loop_hi[0] - loop_lo[0] + 1) / simd<double>::simd_len) - 1);
             // SIMD LOOP
#pragma novector
             for (int x = loop_lo[0]; x <= x_simd_max; x += simd<double>::simd_len)
             {
                 //SIMD compute must be guarded in case compute_t doesn't support it
                 SimdChecker::call_compute<simd<double>>(m_compute,x,y,z);
             }
             // REMAINDER LOOP
#pragma novector
             for (int x = x_simd_max + simd<double>::simd_len; x <= loop_hi[0]; ++x)
             {
                 SimdChecker::call_compute<double>(m_compute,x,y,z);
             }
         }
         else
         {
#warning MK: this is truly horrible code ... the calls are completely unreadable better to have two functions for the inner loop and switch them on or off
             for (int x = loop_lo[0]; x <= loop_hi[0]; ++x)
             {
                 SimdChecker::call_compute<simd<double>>(m_compute,x, y, z);
             }
         }
      }
}

template <class compute_t>
void
FABDriver<compute_t>::execute(const FArrayBox& in, FArrayBox& out)
{
   execute(in,out,out.box());
}

template <class compute_t>
void
FABDriver<compute_t>::execute(const LevelData<FArrayBox>& in, LevelData<FArrayBox>& out, bool fillGhosts)
{
  DataIterator dit0  = in.dataIterator();
  int nbox = dit0.size();
  for(int ibox = 0; ibox < nbox; ++ibox)
  {
    DataIndex di = dit0[ibox];
    const FArrayBox& in_fab = in[di];
    FArrayBox& out_fab = out[di];

    Box out_box;
    if (fillGhosts) out_box = out_fab.box();
    else out_box = in.disjointBoxLayout()[di];

    execute(in_fab,out_fab,out_box);
  }
}

#endif /* FABDRIVER_IMPL_HPP_ */

