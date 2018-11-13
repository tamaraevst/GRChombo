/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

#include "BoxIterator.H"
#include "Copier.H"
#include "FourthOrderInterpStencil.H"
#include "GRLevelData.hpp"
#include "UserVariables.hpp"
#include "DimensionDefinitions.hpp"

/// Class which deals with the boundaries at the edge of the physical domain in
/// cases where they are not periodic. Currently only options are static BCs, 
/// sommerfeld (outgoing radiation) and symmetric. The conditions can differ in 
/// the high and low directions.
class BoundaryConditions
{
  public:
    /// enum for possible boundary states
    enum
    {
        STATIC_BC,
        SOMMERFELD_BC,
        SYMMETRIC_BC
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, CH_SPACEDIM> hi_boundary;
        std::array<int, CH_SPACEDIM> lo_boundary;
    };

    // Member values
    double m_dx;
    int m_num_ghosts;
    params_t m_params;
    ProblemDomain m_domain;
    Box m_domain_box;
    bool is_defined;

    /// Default constructor - need to call define
    BoundaryConditions() { is_defined = false; }

    /// define function sets members and is_defined set to true
    void define(double a_dx, params_t a_params, ProblemDomain a_domain, int a_num_ghosts)
    {
        m_dx = a_dx;
        m_params = a_params;
        m_domain = a_domain;
        m_domain_box = a_domain.domainBox();
        m_num_ghosts = a_num_ghosts;
        is_defined = true;
    }

    /// Fill the rhs boundary values appropriately based on the params set
    void fill_boundary_rhs(const Side::LoHiSide a_side, GRLevelData &a_soln, GRLevelData &a_rhs)
    {
        CH_assert(is_defined);

        // cycle through the directions
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            // only do something if this direction is not periodic
            if (!m_domain.isPeriodic(dir))
            {
                fill_boundary_rhs_dir(a_side, a_soln, a_rhs, dir);
            }
        }
    }

    /// Fill the boundary values appropriately based on the params set
    // in the direction dir
    void fill_boundary_rhs_dir(const Side::LoHiSide a_side, GRLevelData &a_soln,
                                 GRLevelData &a_rhs, const int dir)
    {
        // iterate through the boxes
        for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
        {
            FArrayBox &m_rhs_box = a_rhs[dit()];
            FArrayBox &m_soln_box = a_soln[dit()];
            Box this_box = m_rhs_box.box();
            IntVect offset_lo = - this_box.smallEnd() + m_domain_box.smallEnd();
            IntVect offset_hi = + this_box.bigEnd()   - m_domain_box.bigEnd();

            // reduce box to the intersection of the box and the
            // problem domain ie remove all outer ghost cells
            this_box &= m_domain_box;
            // get the boundary box (may be Empty) and the condition on it
            int boundary_condition = get_boundary_condition(a_side, dir);
            Box boundary_box = get_boundary_box(a_side, dir, offset_lo, offset_hi,
                                                 this_box);

            // now we have the appropriate box, fill it!
            BoxIterator bit(boundary_box);
            for (bit.begin(); bit.ok(); ++bit)
            {
                IntVect iv = bit();
                switch (boundary_condition)
                {
                     case STATIC_BC:
                     {
                         for (int icomp = 0; icomp < NUM_VARS; icomp++)
                         {
                             m_rhs_box(iv, icomp) = 0.0;
                         }
                         break;
                     }
                     case SOMMERFELD_BC:
                     {
                         RealVect loc(iv + 0.5 * RealVect::Unit);
                         loc *= m_dx;
                         double radius_squared = 0.0;
                         FOR1(i) {radius_squared += loc[i]*loc[i];}
                         for (int icomp = 0; icomp < NUM_VARS; icomp++)
                         {
                             m_rhs_box(iv, icomp) = 0.0;
                         }
                         break;
                     }
                     case SYMMETRIC_BC:
                     {
                         IntVect iv_copy = iv;
                         /// where to copy the data from - mirror image in domain
                         if (a_side == Side::Lo)
                         {
                             iv_copy[dir] = -iv[dir] - 1;
                         }
                         else
                         {
                             iv_copy[dir] = 2*m_domain_box.bigEnd(dir) - iv[dir] + 1;
                         }

                         // replace value at iv with value at iv_copy
                         for (int icomp = 0; icomp < NUM_VARS; icomp++)
                         {
                             m_rhs_box(iv, icomp) = m_rhs_box(iv_copy, icomp);
                         }
                         break;
                     }
                     default:
                         MayDay::Error("BoundaryCondition::Supplied boundary not supported.");
                } // end switch
            }     // end iterate over box
        }         // end iterate over boxes
    }

    /// Copy the boundary values from src to dest
    void copy_boundary_cells(const Side::LoHiSide a_side, const GRLevelData &a_src,
                           GRLevelData &a_dest)
    {
        CH_assert(is_defined);
        if(a_src.boxLayout() == a_dest.boxLayout())
        {
            // cycle through the directions
            for (int dir = 0; dir < CH_SPACEDIM; dir++)
            {
                // only do something if this direction is not periodic
                if (!m_domain.isPeriodic(dir))
                {
                    // iterate through the boxes
                    for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
                    {
                        FArrayBox &m_dest_box = a_dest[dit()];
                        Box this_box = m_dest_box.box();
                        IntVect offset_lo = - this_box.smallEnd() + m_domain_box.smallEnd();
                        IntVect offset_hi = + this_box.bigEnd()   - m_domain_box.bigEnd();

                        // reduce box to the intersection of the box and the
                        // problem domain ie remove all outer ghost cells
                        this_box &= m_domain_box;

                        // get the boundary box (may be Empty)
                        Box boundary_box = get_boundary_box(a_side, dir, offset_lo, offset_hi,
                                                 this_box);
            
                        BoxIterator bit(boundary_box);
                        for (bit.begin(); bit.ok(); ++bit)
                        {
                            IntVect iv = bit();
                            for (int icomp = 0; icomp < NUM_VARS; icomp++)
                            {
                                m_dest_box(iv, icomp) = a_src[dit()](iv, icomp);
                            }
                        }         // end iterate over box
                    }             // end iterate over boxes
                }                 // end if(not periodic)
            }                     // end iterate over spacedims
        }                         // end test for same box layout
    }

    /// enforce symmetric boundary conditions, e.g. after interpolation
    void enforce_symmetric_boundaries(const Side::LoHiSide a_side, GRLevelData &a_state)
    {
        CH_assert(is_defined);

        // cycle through the directions
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            // only do something if this direction is not periodic and symmetric
            if (!m_domain.isPeriodic(dir))
            {
                int boundary_condition = get_boundary_condition(a_side, dir);

                // as a bit of a hack, just use the rhs update, since it is the same copying of cells
                if (boundary_condition == SYMMETRIC_BC) 
                {
                    fill_boundary_rhs_dir(a_side, a_state, a_state, dir);
                }
            }
        }
    }

    /// Fill the fine boundary values in a_state 
    void interp_boundaries(GRLevelData &a_fine_state, 
                           GRLevelData &a_coarse_state, const Side::LoHiSide a_side)
    {
        CH_assert(is_defined);

        // cycle through the directions
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            // only do something if this direction is not periodic
            if (!m_domain.isPeriodic(dir))
            {
                // Ref ratio is always two
                int ref_ratio = 2;

                // create a coarsened fine layout and copy the coarse data onto it
                DisjointBoxLayout coarsened_layout;
                coarsen(coarsened_layout, a_fine_state.disjointBoxLayout(), ref_ratio*IntVect::Unit);
                GRLevelData coarsened_fine;
                coarsened_fine.define(coarsened_layout, NUM_VARS, m_num_ghosts*IntVect::Unit);
                Box coarse_domain_box = coarsen(m_domain_box, ref_ratio);

                // trick the copyTo into thinking the boundary cells are within the domain by growing the domain
                Box grown_domain_box = coarse_domain_box;
                grown_domain_box.grow(m_num_ghosts*IntVect::Unit);
                Copier boundary_copier;
                boundary_copier.ghostDefine(a_coarse_state.disjointBoxLayout(), coarsened_fine.disjointBoxLayout(),
                          grown_domain_box, m_num_ghosts*IntVect::Unit, m_num_ghosts*IntVect::Unit);
                a_coarse_state.copyTo(a_coarse_state.interval(), coarsened_fine, 
                                      coarsened_fine.interval(), boundary_copier);

                // iterate through the boxes, interpolating onto the fine cells
                for (DataIterator dit = coarsened_layout.dataIterator(); dit.ok(); ++dit)
                {
                    FArrayBox &m_fine_box = a_fine_state[dit()];
                    FArrayBox &m_coarse_box = coarsened_fine[dit()];
                    Box this_box = m_coarse_box.box();
                    Box fine_box = m_fine_box.box();
                    IntVect offset_lo = - this_box.smallEnd() + coarse_domain_box.smallEnd();
                    IntVect offset_hi = + this_box.bigEnd()   - coarse_domain_box.bigEnd();
        
                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= coarse_domain_box;

                    // get the boundary box - remove one cell as we only want 2 coarse cells
                    // filled in each direction to fill the 3 fine cells above
                    Box boundary_box = get_boundary_box(a_side, dir, offset_lo, offset_hi,
                                                 this_box, 1);

                    //define standard stencil for interp where not near boundaries in other dirs
                    IntVect default_offset = IntVect::Zero;
                    if (a_side == Side::Hi) {default_offset[dir] = +2;}
                    if (a_side == Side::Lo) {default_offset[dir] = -2;}
                    FourthOrderInterpStencil default_stencil(default_offset, ref_ratio);

                    // now interp the box
                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        IntVect iv = bit();
                        IntVect lo_local_offset = iv - m_coarse_box.smallEnd();
                        IntVect hi_local_offset = m_coarse_box.bigEnd() - iv;

                        // bit of work to get the right stencils for near the edges of the box
                        bool near_boundary = false;
                        IntVect local_boundary_offset = IntVect::Zero;
                        for (int dir2=0; dir2 < CH_SPACEDIM; dir2++)
                        {
                            if(dir2 == dir)
                            {
                                local_boundary_offset[dir2] = default_offset[dir2];
                            }
                            else if((dir2 !=dir) && (lo_local_offset[dir2] > 1) 
                                        && (hi_local_offset[dir2] > 1))
                            {
                                local_boundary_offset[dir2] = 0;
                            }
                            else if((dir2 != dir) && (lo_local_offset[dir2] == 1))
                            {
                                local_boundary_offset[dir2] = -2;
                                near_boundary = true;
                            }
                            else if((dir2 !=dir) && (hi_local_offset[dir2] == 1))
                            {
                                local_boundary_offset[dir2] = +2;
                                near_boundary = true;
                            }
                            else
                            {
                                MayDay::Error("BoundaryConditions::define bad boxes");
                            }
                        }

                        // if not near the boundary use the default stencil, otherwise
                        // use the one calculated locally
                        if (!near_boundary)
                        {
                            default_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                        }
                        else
                        {
                            FourthOrderInterpStencil local_stencil(local_boundary_offset, ref_ratio);
                            local_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                        }
                    }  // end loop box
                }      // end loop boxes
            }          // end if isPeriodic
        }              // end loop dir
    }

    /// Get the boundary condition for dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int dir)
    {
        int boundary_condition = 0;
        if (a_side == Side::Lo)
        {
            boundary_condition = m_params.lo_boundary[dir];
        }
        else
        {
            boundary_condition = m_params.hi_boundary[dir];
        }
        return boundary_condition;
    }

    /// get the boundary box to fill if we are at a boundary
    Box get_boundary_box(const Side::LoHiSide a_side, const int dir, const IntVect &offset_lo,
                   const IntVect &offset_hi, Box &this_ghostless_box, int shrink_for_coarse = 0)
    {
        // default constructor gives empty box
        Box boundary_box;

        // check if we are over the edges of the domain - are we a boundary box?
        // if so create the box of the cells we want to fill
        if (((a_side == Side::Hi) && (offset_hi[dir] > 0)) || 
            ((a_side == Side::Lo) && (offset_lo[dir] > 0)))
        {
            // Get just the boundary box to iterate over, m_num_ghosts ghost
            // cells unless we are filling the coarse cells in the interp case
            // where we want to fill only two coarse ghost cells (to cover 3 fine ones)
            if (a_side == Side::Lo)
            {
                boundary_box = adjCellLo(this_ghostless_box, dir, m_num_ghosts - shrink_for_coarse);
            }
            else
            {
                boundary_box = adjCellHi(this_ghostless_box, dir, m_num_ghosts - shrink_for_coarse);
            }

            // adjust for any offsets - catches the corners etc
            // but only want to fill them once, so x fills y and z, y fills z etc
            for (int dir2 = 0 ; dir2 < CH_SPACEDIM ; dir2 ++)
            {
                if (offset_lo[dir2] > 0) // this direction is a low end boundary
                {
                    if ((dir2 > dir) && !(m_domain.isPeriodic(dir2)))
                    {
                        //grow it to fill the corners
                        boundary_box.growLo(dir2, m_num_ghosts - shrink_for_coarse);
                    }
                }
                else // cut off end ghost cell
                {
                    if (dir2 != dir) {boundary_box.growLo(dir2, - shrink_for_coarse);}
                }

                if (offset_hi[dir2] > 0) //this direction is a high end boundary
                {
                    if ((dir2 > dir) && !(m_domain.isPeriodic(dir2)))
                    {
                        //grow it to fill the corners
                        boundary_box.growHi(dir2, m_num_ghosts - shrink_for_coarse);
                    }
                }
                else // cut off end ghost cell
                {
                    if (dir2 != dir) {boundary_box.growHi(dir2, - shrink_for_coarse);}
                }
            }
        }
        return boundary_box;
    }
};

#endif /* BOUNDARYCONDITIONS_HPP_ */
