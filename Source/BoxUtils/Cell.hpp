#ifndef CELL_HPP_
#define CELL_HPP_

#include "AlwaysInline.hpp"
#include "BoxPointers.hpp"
#include "CellIndex.hpp"
#include "GRInterval.hpp"
#include "IntVect.H"
#include "tensor.hpp"

/// Encapsulates information about the position of a cell
/** It contains the position of the cell on the Chombo grid and the index of the
 * flattened Chombo array where the data is stored (both for the input and the
 * output data of the FABDriver). The main use of Cell is to hand information on
 * the current cells from the FABDriver to the compute classes.
 */
template <class data_t> class Cell
{
  protected:
    const IntVect
        m_integer_coords; //!< Integer coordinates of the cell in the grid

    /// Index in the flattened Chombo array where the input data for this cell
    /// resides
    const CellIndexIn m_in_index;
    /// Index in the flattened Chombo array where the output data for this cell
    /// should be written to
    const CellIndexOut m_out_index;

    /// Contains pointers to the low and high ends of the current box and the
    /// strides
    const BoxPointers &m_box_pointers;

  public:
    Cell(const IntVect integer_coords, const BoxPointers &box_pointers)
        : m_integer_coords(integer_coords),
          m_in_index(box_pointers.get_in_index(integer_coords)),
          m_out_index(box_pointers.get_out_index(integer_coords)),
          m_box_pointers(box_pointers)
    {
    }

    // Don't accept rvalues in constructor since we only store reference to
    // box_pointers
    Cell(const IntVect integer_coords, BoxPointers &&box_pointers) = delete;

    /// Allows implicit conversion from Cell to CellIndexIn.
    //! As a result of this definition one may pass a Cell to a function
    //! expecting a CellIndexIn
    ALWAYS_INLINE
    operator CellIndexIn() const { return m_in_index; }

    /// Allows implicit conversion from Cell to CellIndexOut.
    //! As a result of this definition one may pass a Cell to a function
    //! expecting a CellIndexOut
    ALWAYS_INLINE
    operator CellIndexOut() const { return m_out_index; }

    /// Allows implicit conversion from Cell to IntVect
    ALWAYS_INLINE
    operator IntVect() const { return m_integer_coords; }

    ALWAYS_INLINE
    bool operator==(const IntVect iv) { return (m_integer_coords == iv); }

    /// Returns the integer coordinates of this Cell in the form of an IntVect
    ALWAYS_INLINE
    IntVect get_int_vect() const { return m_integer_coords; }

    /// Returns Index in the flattened Chombo array where the input data for
    /// this cell resides
    ALWAYS_INLINE
    CellIndexIn get_in_index() const { return m_in_index; }

    /// Returns Index in the flattened Chombo array where the output data for
    /// this cell should be written to
    ALWAYS_INLINE
    CellIndexOut get_out_index() const { return m_out_index; }

    D_TERM( // Chombo's dimension independent macro. We only want iy() and iz()
            // if we run in enough dimensions
        ALWAYS_INLINE int ix() const { return m_integer_coords[0]; },
        ALWAYS_INLINE int iy() const { return m_integer_coords[1]; },
        ALWAYS_INLINE int iz() const { return m_integer_coords[2]; })

    /// Returns the box pointers
    ALWAYS_INLINE
    const BoxPointers &get_box_pointers() const { return m_box_pointers; }

    ALWAYS_INLINE
    data_t load_vars(int icomp) const;

    ALWAYS_INLINE
    void load_vars(data_t &out, int icomp) const;

    void load_vars(data_t (&out)[NUM_VARS]) const;

    template <template <typename> class vars_t>
    void load_vars(vars_t<data_t> &vars) const;

    template <template <typename> class vars_t> auto load_vars() const;

    void store_vars(const data_t &value, const int icomp) const;

    template <int start_var, int end_var>
    void store_vars(
        const tensor<1, data_t, GRInterval<start_var, end_var>::size()> &values,
        const GRInterval<start_var, end_var> interval) const;

    void store_vars(const std::array<data_t, NUM_VARS> &values) const;

    template <template <typename> class vars_t>
    void store_vars(vars_t<data_t> &vars) const;
};

#include "Cell.impl.hpp"

#endif /* CELL_HPP_ */
