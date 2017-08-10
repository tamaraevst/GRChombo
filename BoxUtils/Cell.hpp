#ifndef CELL_HPP_
#define CELL_HPP_

#include "always_inline.hpp"
#include "CellIndex.hpp"
#include "BoxPointers.hpp"
#include "IntVect.H"
#include "VarsBase.hpp"

///Encapsulates information about the position of a cell
/** It contains the position of the cell on the Chombo grid and the index of the flattened Chombo array where the data
  * is stored (both for the input and the output data of the FABDriver).
  * The main use of Cell is to hand information on the current cells from the FABDriver to the compute classes.
  */
template <class data_t>
class Cell
{
protected:
    const IntVect m_integer_coords; //!< Integer coordinates of the cell in the grid

    /// Index in the flattened Chombo array where the input data for this cell resides
    const CellIndexIn m_in_index;
    ///Index in the flattened Chombo array where the output data for this cell should be written to
    const CellIndexOut m_out_index;

    ///Contains pointers to the low and high ends of the current box and the strides
    const BoxPointers& m_box_pointers;

public:
    Cell(const IntVect integer_coords, const BoxPointers& box_pointers) :
        m_integer_coords (integer_coords),
        m_in_index( box_pointers.get_in_index(integer_coords) ),
        m_out_index( box_pointers.get_out_index(integer_coords) ),
        m_box_pointers (box_pointers)
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

    ALWAYS_INLINE
    bool operator ==(const IntVect iv)
    {
        return (m_integer_coords == iv);
    }

    ///Returns the integer coordinates of this Cell in the form of an IntVect
    ALWAYS_INLINE
    IntVect get_int_vect() const
    {
        return m_integer_coords;
    }

    ///Returns Index in the flattened Chombo array where the input data for this cell resides
    ALWAYS_INLINE
    CellIndexIn get_in_index() const
    {
        return m_in_index;
    }

    ///Returns Index in the flattened Chombo array where the output data for this cell should be written to
    ALWAYS_INLINE
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

    ///Returns the box pointers
    ALWAYS_INLINE
    BoxPointers get_box_pointers() const
    {
        return m_box_pointers;
    }

    ALWAYS_INLINE
    data_t local_vars(int icomp) const;

    ALWAYS_INLINE
    void local_vars(data_t& out, int icomp) const;

    void local_vars(data_t (&out)[c_NUM]) const;

    template<template<typename> class vars_t>
    void local_vars(vars_t<data_t>& vars) const;

    template<template<typename> class vars_t>
    vars_t<data_t> local_vars() const;

    void store_vars(const data_t& value, const int icomp) const;

    template <int num_comp>
    void store_vars(const tensor<1, data_t, num_comp>& values, const int start_comp) const;

    void store_vars(const std::array<data_t, c_NUM>& values) const;

    template<template<typename> class vars_t>
    void store_vars(vars_t<data_t>& vars) const;
};

#include "Cell.impl.hpp"

#endif /* CELL_HPP_ */
