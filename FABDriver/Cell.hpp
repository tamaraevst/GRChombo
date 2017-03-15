#ifndef CELL_HPP_
#define CELL_HPP_

#include "CellIndex.hpp"

class Cell
{
public:

    const int m_ix;
    const int m_iy;
    const int m_iz;

    const CellIndexIn m_in_index;
    const CellIndexOut m_out_index;

    Cell(const int ix, const int iy, const int iz, const int in_lo[CH_SPACEDIM], const int out_lo[CH_SPACEDIM],
                 const int  in_stride[CH_SPACEDIM], const int out_stride[CH_SPACEDIM]):
        m_ix(ix), m_iy(iy), m_iz(iz),
        m_in_index(in_stride[2]*(iz-in_lo[2]) + in_stride[1]*(iy-in_lo[1]) + (ix-in_lo[0])),
        m_out_index(out_stride[2]*(iz-out_lo[2]) + out_stride[1]*(iy-out_lo[1]) + (ix-out_lo[0]))
    {}

    operator CellIndexIn () const
    {
        return m_in_index;
    }

    operator CellIndexOut () const
    {
        return m_out_index;
    }
};

#endif /* CELL_HPP_ */
