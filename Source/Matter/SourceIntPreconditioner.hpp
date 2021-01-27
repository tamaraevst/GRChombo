/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SOURCEINTPRECONDITIONER_HPP
#define SOURCEINTPRECONDITIONER_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "MatterCCZ4.hpp"
#include "simd.hpp"

//! Sets a field variable to zero outside a set radius
//! constructor on the grid
template <class matter_t> class SourceIntPreconditioner
{
  public:
    template <class data_t>
    using Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;

    //! Constructor
    SourceIntPreconditioner(const matter_t &a_matter, const double dx, const double a_L,
             const std::array<double,CH_SPACEDIM> a_centre,
             const int a_c_var = -1, const double a_radius = 0.)
             : m_matter(a_matter), m_deriv(dx),  m_dx(dx), m_L(a_L),
             m_centre(a_centre), m_c_var(a_c_var), m_radius(a_radius) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx,m_centre);
        const auto vars = current_cell.template load_vars<Vars>();

        data_t r = coords.get_radius();
        data_t r_max = m_radius;
        data_t heaviside = 0.5*(1.- abs(r-r_max)/(r-r_max));
        auto original_val = current_cell.load_vars(m_c_var);
        data_t output_val = original_val*heaviside;
        current_cell.store_vars(output_val, m_c_var);

    }

  protected:

    const double m_L, m_dx, m_radius; // simulation boxsize and cellsize
    const matter_t &m_matter;
    FourthOrderDerivatives m_deriv;
    const int m_c_var;      // var enum for the angmom source
    const std::array<double,CH_SPACEDIM> m_centre; // centre of mom flux calc
};

//#include "SourceIntPreconditioner.impl.hpp"

#endif /* SOURCEINTPRECONDITIONER_HPP */
