#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

#include "simd.hpp"

template <class data_t>
class Coordinates
{
public:
    data_t x; //We vectorise over x so we must allow x to be a vector
    double y;
    double z;

    Coordinates(IntVect integer_coords, double dx)
    {
        compute_coord<data_t>(x, integer_coords[0], dx);

        //The below code allows for 2D Cartoon reduction:
#if IDX_SPACEDIM == CH_SPACEDIM && CH_SPACEDIM == 3
        compute_coord<double>(y, integer_coords[1], dx);
        compute_coord<double>(z, integer_coords[2], dx);
#elif IDX_SPACEDIM == CH_SPACEDIM + 1 && CH_SPACEDIM == 2
        y = 0;
        compute_coord<double>(z, integer_coords[1], dx);
#else
#ifdef CH_SPACEDIM
#error compute_coord has not got your dimension combination implemented.
#endif
#endif
    }

    template <typename t>
    ALWAYS_INLINE
    static
    void
    compute_coord(t& out, int position, double dx)
    {
        out = (position+0.5)*dx;
    }

    //MK: I passed 'out' as argument because overloading by return type doesn't work
    ALWAYS_INLINE
    static
    typename std::enable_if<(simd_traits<double>::simd_len > 1), void >::type
    compute_coord(simd<double>& out, int position, double dx)
    {
        double out_arr[simd_traits<double>::simd_len];
        for (int i = 0; i < simd_traits<double>::simd_len; ++i)
        {
            out_arr[i] = (position+i+0.5)*dx;
        }
        out = simd<double>::load(out_arr);
    }
};
#endif /* COORDINATES_HPP_ */
