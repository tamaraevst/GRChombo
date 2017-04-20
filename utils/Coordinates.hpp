#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

#include "simd.hpp"
#include "Cell.hpp"

template <class data_t>
class Coordinates
{
public:
    data_t x; //We vectorise over x so we must allow x to be a vector
    double y;
    double z;

    Coordinates(IntVect integer_coords, double dx)
    {
        compute_coord(x, integer_coords[0], dx);

        //The below code allows for 2D Cartoon reduction:
#if IDX_SPACEDIM == CH_SPACEDIM && CH_SPACEDIM == 3
        compute_coord(y, integer_coords[1], dx);
        compute_coord(z, integer_coords[2], dx);
#elif IDX_SPACEDIM == CH_SPACEDIM + 1 && CH_SPACEDIM == 2
        y = 0;
        compute_coord(z, integer_coords[1], dx);
#else
#ifdef CH_SPACEDIM
#error compute_coord has not got your dimension combination implemented.
#endif
#endif
    }

    ALWAYS_INLINE
    static
    void
    compute_coord(double& out, int position, double dx)
    {
        out = (position+0.5)*dx;
    }

    //MK: I passed 'out' as argument because overloading by return type doesn't work
    ALWAYS_INLINE
    static
    void
    compute_coord(simd<double>& out, int position, double dx)
    {
        double out_arr[simd_traits<double>::simd_len];
        for (int i = 0; i < simd_traits<double>::simd_len; ++i)
        {
            out_arr[i] = (position+i+0.5)*dx;
        }
        out = simd<double>::load(out_arr);
    }

    data_t
    get_radius(IntVect center)
    {
        //Note that this is not currently dimension independent
	data_t r = sqrt(  pow(x - center[0],2)
                        + pow(y - center[1],2)
                        + pow(z - center[2],2));

        double minimum_r = 1e-6;
        auto r_is_too_small = simd_compare_lt(r, minimum_r);
        return simd_conditional(r_is_too_small, minimum_r, r);
    }

    static
    data_t
    get_radius(IntVect integer_coords, double dx, IntVect center)
    {
        data_t xx;
        double yy;
        double zz;

        //Note that this is not currently dimension independent
        compute_coord(xx, integer_coords[0], dx);
        compute_coord(yy, integer_coords[1], dx);
        compute_coord(zz, integer_coords[2], dx);
 
	data_t r = sqrt(  pow(xx - center[0],2)
                        + pow(yy - center[1],2)
                        + pow(zz - center[2],2));

        double minimum_r = 1e-6;
        auto r_is_too_small = simd_compare_lt(r, minimum_r);
        return simd_conditional(r_is_too_small, minimum_r, r);
    }

};
#endif /* COORDINATES_HPP_ */
