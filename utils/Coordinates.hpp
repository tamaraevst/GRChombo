#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

template <class data_t>
class Coordinates
{
public:
    data_t x; //We vectorise over x so we must allow x to be a vector
    double y;
    double z;

    Coordinates(int ix, int iy, int iz, double dx)
    {
        compute_coord<data_t>(x, ix, dx);

        //The below code allows for 2D Cartoon reduction:
#if IDX_SPACEDIM == CH_SPACEDIM
        compute_coord<double>(y, iy, dx);
#elif IDX_SPACEDIM == CH_SPACEDIM + 1
        y = 0;
#else
#ifdef CH_SPACEDIM
#error compute_coord has not got your dimension combination implemented.
#endif
#endif

        compute_coord<double>(z, iz, dx);
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
