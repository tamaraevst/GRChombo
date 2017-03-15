#ifndef CELL_HPP_
#define CELL_HPP_

#include "CellIndex.hpp"

///Encapsulates information about the position of a cell
/** It contains the position of the cell on the Chombo grid and the index of the flattened Chombo array where the data
  * is stored (both for the input and the output data of the FABDriver).
  * The main use of Cell is to hand information on the current cells from the FABDriver to the compute classes.
  */
class Cell
{
public:
    const int m_ix; //!< Integer x coordinate of the cell in the grid
    const int m_iy; //!< Integer y coordinate of the cell in the grid
    const int m_iz; //!< Integer z coordinate of the cell in the grid

    /// Index in the flattened Chombo array where the input data for this cell resides
    const CellIndexIn m_in_index;
    ///Index in the flattened Chombo array where the output data for this cell should be written to
    const CellIndexOut m_out_index;

    Cell(const int ix, const int iy, const int iz, const int in_lo[CH_SPACEDIM], const int out_lo[CH_SPACEDIM],
                 const int in_stride[CH_SPACEDIM], const int out_stride[CH_SPACEDIM]):
        m_ix(ix), m_iy(iy), m_iz(iz),
        m_in_index(in_stride[2]*(iz-in_lo[2]) + in_stride[1]*(iy-in_lo[1]) + (ix-in_lo[0])),
        m_out_index(out_stride[2]*(iz-out_lo[2]) + out_stride[1]*(iy-out_lo[1]) + (ix-out_lo[0]))
    {}

    ///Allows implicit conversion from Cell to CellIndexIn.
    //!As a result of this definition one may pass a Cell to a function expecting a CellIndexIn
    ALWAYS_INLINE
    operator CellIndexIn () const
    {
        return m_in_index;
    }

    ///Allows implicit conversion from Cell to CellIndexOut.
    //!As a result of this definition one may pass a Cell to a function expecting a CellIndexOut
    ALWAYS_INLINE
    operator CellIndexOut () const
    {
        return m_out_index;
    }
};

#endif /* CELL_HPP_ */
