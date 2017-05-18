// Last edited K Clough 17.05.17

#if !defined(SCALARFIELDLEVEL_HPP_)
#error "This file should only be included through ScalarFieldLevel.hpp"
#endif

#ifndef SCALARFIELDLEVEL_IMPL_HPP_
#define SCALARFIELDLEVEL_IMPL_HPP_

//General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"

//For RHS update
#include "CCZ4Matter.hpp"

//For constraints calculation
#include "ConstraintsMatter.hpp"

//For tag cells
#include "ComputeModGrad.hpp"

//Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "ScalarBubble.hpp"
#include "RelaxationChi.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    //Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS, no_simd_support());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity) pout () << "ScalarFieldLevel::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc, then
    //initial conditions for scalar field - here a bubble
    BoxLoops::loop(make_compute_pack(SetValue(0.0), ScalarBubble(m_p.initial_params, m_dx)), m_state_new, m_state_new, FILL_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(ConstraintsMatter<ScalarFieldWithPotential>(scalar_field, m_dx,
                    m_p.G_Newton), m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{

    //Relaxation function for chi - this will eventually be done separately with hdf5 as input
    if (m_time < m_p.relaxtime)
    {
        //Calculate chi relaxation right hand side
        //Note this assumes conformal chi and Mom constraint trivially satisfied
        //No evolution in other variables, which are assumed to satisfy constraints per initial conditions
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(make_compute_pack(RelaxationChi<ScalarFieldWithPotential>(scalar_field, m_dx, m_p.relaxspeed,
                                   m_p.G_Newton), SetValue(0.0, Interval(c_h11, c_Mom3))), a_soln, a_rhs, SKIP_GHOST_CELLS);
    }
    else
    {

        //Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(make_compute_pack(EnforceTfA(), PositiveChiAndAlpha()), a_soln, a_soln, FILL_GHOST_CELLS);

        //Calculate CCZ4Matter right hand side with matter_t = ScalarField
        //We don't want undefined values floating around in the constraints so zero these
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        BoxLoops::loop(make_compute_pack(CCZ4Matter<ScalarFieldWithPotential>(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                        m_p.formulation, m_p.G_Newton), SetValue(0.0, Interval(c_Ham, c_Mom3))), a_soln, a_rhs, SKIP_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce trace free A_ij
    BoxLoops::loop(EnforceTfA(), a_soln, a_soln, FILL_GHOST_CELLS);
}

// Specify if you want any plot files to be written, with which vars
void ScalarFieldLevel::specificWritePlotHeader(std::vector<int> &plot_states) const
{
    plot_states = {c_chi, c_K};
}

// Tell Chombo when to mark cells for regridding
// override virtual tagcells function for tagging based on K and phi gradients
// TODO KClough: It is a bit ugly to have all this here, may want to create a hook
// in the main routine and then just put the conditional here - but tagging is very
// problem specific, so it is hard to have a "general" routine
void ScalarFieldLevel::tagCells (IntVectSet& a_tags)
{
    CH_TIME("GRAMRLevel::tagCells");
    if ( m_verbosity ) pout () << "GRAMRLevel::tagCells " << m_level << endl;

    fillAllGhosts(); //We need filled ghost cells to calculate gradients etc

    // Create tags based on undivided gradient of phi and K
    IntVectSet local_tags;

    const DisjointBoxLayout& level_domain = m_state_new.disjointBoxLayout();
    DataIterator dit0 = level_domain.dataIterator();
    int nbox = dit0.size();
    for(int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        const Box& b = level_domain[di];
        const FArrayBox& state_fab = m_state_new[di];

        //mod gradient
        FArrayBox mod_grad_fab(b,c_NUM);
        BoxLoops::loop(ComputeModGrad(m_dx), state_fab, mod_grad_fab);

        const IntVect& smallEnd = b.smallEnd();
        const IntVect& bigEnd = b.bigEnd();

        const int xmin = smallEnd[0];
        const int ymin = smallEnd[1];
        const int zmin = smallEnd[2];

        const int xmax = bigEnd[0];
        const int ymax = bigEnd[1];
        const int zmax = bigEnd[2];

#pragma omp parallel for collapse(3) schedule(static) default(shared)
        for (int iz = zmin; iz <= zmax; ++iz)
        for (int iy = ymin; iy <= ymax; ++iy)
        for (int ix = xmin; ix <= xmax; ++ix)
        {
            IntVect iv(ix,iy,iz);

            //At the moment only base on gradient K and phi
            //NB the mod_grad_fab is the true gradient, but we want the change
            //across the cell, so multiply by dx (otherwise you regrid forever...)
            if ((mod_grad_fab(iv,c_K)* m_dx >= m_p.regrid_threshold_K)
            || (mod_grad_fab(iv,c_phi) * m_dx >= m_p.regrid_threshold_phi))
            {
                // local_tags |= is not thread safe.
#pragma omp critical
                {
                    local_tags |= iv;
                }
            }
        }
    }

    local_tags.grow(m_p.tag_buffer_size);

    // Need to do this in two steps unless a IntVectSet::operator &=
    // (ProblemDomain) operator is defined
    Box local_tags_box = local_tags.minBox();
    local_tags_box &= m_problem_domain;
    local_tags &= local_tags_box;

    a_tags = local_tags;
}

#endif
