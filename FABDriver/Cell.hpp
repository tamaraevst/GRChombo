#ifndef CELL_HPP_
#define CELL_HPP_

#include "CellIndex.hpp"
#include "IntVect.H"

///Encapsulates information about the position of a cell
/** It contains the position of the cell on the Chombo grid and the index of the flattened Chombo array where the data
  * is stored (both for the input and the output data of the FABDriver).
  * The main use of Cell is to hand information on the current cells from the FABDriver to the compute classes.
  */
class Cell
{
protected:
    const IntVect m_integer_coords; //!< Integer coordinates of the cell in the grid

    /// Index in the flattened Chombo array where the input data for this cell resides
    const CellIndexIn m_in_index;
    ///Index in the flattened Chombo array where the output data for this cell should be written to
    const CellIndexOut m_out_index;

public:
    Cell(const IntVect integer_coords, const int in_lo[CH_SPACEDIM], const int out_lo[CH_SPACEDIM],
                 const int in_stride[CH_SPACEDIM], const int out_stride[CH_SPACEDIM]):
        m_integer_coords (integer_coords),
        m_in_index( in_stride[2]*(integer_coords[2]-in_lo[2]) +
                    in_stride[1]*(integer_coords[1]-in_lo[1]) + (integer_coords[0]-in_lo[0]) ),
        m_out_index( out_stride[2]*(integer_coords[2]-out_lo[2]) +
                     out_stride[1]*(integer_coords[1]-out_lo[1]) + (integer_coords[0]-out_lo[0]) )
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

    ///Allows implicit conversion from Cell to IntVect
    ALWAYS_INLINE
    operator IntVect () const
    {
        return m_integer_coords;
    }

    bool operator ==(const IntVect iv)
    {
        return (m_integer_coords == iv);
    }

    ///Returns the integer coordinates of this Cell in the form of an IntVect
    IntVect get_int_vect() const
    {
        return m_integer_coords;
    }

    ///Returns Index in the flattened Chombo array where the input data for this cell resides
    CellIndexIn get_in_index() const
    {
        return m_in_index;
    }

    ///Returns Index in the flattened Chombo array where the output data for this cell should be written to
    CellIndexOut get_out_index() const
    {
        return m_out_index;
    }

    D_TERM( //Chombo's dimension independent macro. We only want iy() and iz() if we run in enough dimensions
    ALWAYS_INLINE
    int ix() const { return m_integer_coords[0]; },
    ALWAYS_INLINE
    int iy() const { return m_integer_coords[1]; },
    ALWAYS_INLINE
    int iz() const { return m_integer_coords[2]; }
    )
};

#endif /* CELL_HPP_ */
