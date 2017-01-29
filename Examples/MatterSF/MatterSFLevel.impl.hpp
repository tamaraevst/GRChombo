// Last edited K Clough 16.10.16

#if !defined(MATTERSFLEVEL_HPP_)
#error "This file should only be included through MatterSFLevel.hpp"
#endif

#ifndef MATTERSFLEVEL_IMPL_HPP_
#define MATTERSFLEVEL_IMPL_HPP_

//General includes common to most GR problems
#include "MatterSFLevel.hpp"
#include "FABDriver.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"

//For RHS update
//#include "CCZ4.hpp"
#include "CCZ4Matter.hpp"

//For constraints calculation
#include "ConstraintsMatter.hpp"

//For tag cells
#include "ComputeModGrad.hpp"

//Problem specific includes
#include "BubbleSF.hpp"
#include "SFMatter.hpp"
#include "RelaxationChi.hpp"

// Things to do at each advance step, after the RK4 is calculated
void MatterSFLevel::specificAdvance()
{
    //Enforce the trace free A_ij condition
    FABDriver<EnforceTfA>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) FABDriver<NanCheck>().execute(m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void MatterSFLevel::initialData()
{
    CH_TIME("MatterSFLevel::initialData");
    if (m_verbosity) pout () << "MatterSFLevel::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc
    m_state_new.setVal(0.);

    //Initial conditions for scalar field - here a bubble
    FABDriver<BubbleSF>(m_p.matter_params, m_dx).execute(m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());

}

// Things to do before outputting a checkpoint file
void MatterSFLevel::preCheckpointLevel()
{
    fillAllGhosts();
    FABDriver<ConstraintsMatter<SFMatter> >(m_p.matter_params, m_dx, m_p.G_Newton).execute(m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void MatterSFLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{

    //Relaxation function for chi - this will eventually be done separately with hdf5 as input
    if (m_time < m_p.relaxtime) {
      //Calculate chi relaxation right hand side
      //Note this assumes conformal chi and Mom constraint trivially satisfied
      FABDriver<RelaxationChi<SFMatter> >(m_p.matter_params, m_dx, m_p.relaxspeed, m_p.G_Newton).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

       //No evolution in other variables, which are assumed to satisfy constraints per initial conditions
      a_rhs.setVal(0., Interval(c_h11,c_Mom3));
    } else {

    	FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

   	 	//Enforce positive chi and alpha
    	FABDriver<PositiveChiAndAlpha>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

    	//Calculate CCZ4 right hand side with SF matter
    	FABDriver<CCZ4Matter<SFMatter> >(m_p.ccz4Params, m_p.matter_params, m_dx, m_p.sigma, m_p.formulation, m_p.G_Newton).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

    	//We don't want undefined values floating around in the constraints
    	a_rhs.setVal(0., Interval(c_Ham,c_Mom3));

    }
}

// Things to do at ODE update, after soln + rhs
void MatterSFLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free A_ij condition
    FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);
}

// Tell Chombo when to mark cells for regridding
// override virtual tagcells function for tagging based on K and phi gradients
// TODO KClough: It is a bit ugly to have all this here, may want to create a hook
// in the main routine and then just put the conditional here - but tagging is very
// problem specific, so it is hard to have a "general" routine
void MatterSFLevel::tagCells (IntVectSet& a_tags)
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
        FABDriver<ComputeModGrad>(m_dx).execute(state_fab, mod_grad_fab);

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
