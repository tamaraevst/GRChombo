#if !defined(MATTERSFLEVEL_HPP_)
#error "This file should only be included through MatterSFLevel.hpp"
#endif

#ifndef MATTERSFLEVEL_IMPL_HPP_
#define MATTERSFLEVEL_IMPL_HPP_

#include "MatterSFLevel.hpp"
#include "FABDriver.hpp"
#include "EnforceTfA.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "NanCheck.hpp"
#include "ConstraintsMatter.hpp"
#include "CCZ4SFMatter.hpp"
#include "RelaxationChi.hpp"

//Initial data
#include "BubbleSF.hpp"

void MatterSFLevel::specificAdvance()
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(m_state_new, m_state_new, FILL_GHOST_CELLS);

    //Check for nan's
    if (m_p.nan_check) FABDriver<NanCheck>().execute(m_state_new, m_state_new, SKIP_GHOST_CELLS, disable_simd());
}

void MatterSFLevel::initialData()
{
    CH_TIME("MatterSFLevel::initialData");
    if (m_verbosity) pout () << "MatterSFLevel::initialData " << m_level << endl;

    //First set everything to zero ... we don't want undefined values in constraints etc
    m_state_new.setVal(0.);

    //Initial conditions for scalar field - here a bubble
    FABDriver<BubbleSF>(m_p.sfm_params, m_dx).execute(m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());

}

void MatterSFLevel::preCheckpointLevel()
{
    fillAllGhosts();
    FABDriver<ConstraintsMatter>(m_dx).execute(m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

void MatterSFLevel::specificEvalRHS(GRLevelData& a_soln, GRLevelData& a_rhs, const double a_time)
{

    //Relaxation function for chi - this will eventually be done separately with hdf5 as input
    if (m_time < m_p.relaxtime) {

       //Calculate chi relaxation right hand side
       //Note that this assumes conformal chi and that the momentum constraint is trivially satisfied
       FABDriver<RelaxationChi>(m_dx, m_p.relaxspeed).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

       //No evolution in other variables, which are assumed to satisfy constraints per initial conditions
       a_rhs.setVal(0., Interval(c_h11,c_Mom3));
    } 

    //Else do normal CCZ4 evolution
    else {

    	FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

   	 	//Enforce positive chi and alpha
    	FABDriver<PositiveChiAndAlpha>().execute(a_soln, a_soln, FILL_GHOST_CELLS);

    	//Calculate CCZ4 right hand side with SF matter
    	FABDriver<CCZ4SFMatter>(m_p.ccz4Params, m_dx, m_p.sigma).execute(a_soln, a_rhs, SKIP_GHOST_CELLS);

    	//We don't want undefined values floating around in the constraints
    	a_rhs.setVal(0., Interval(c_Ham,c_Mom3));

    }
}

void MatterSFLevel::specificUpdateODE(GRLevelData& a_soln, const GRLevelData& a_rhs, Real a_dt)
{
    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(a_soln, a_soln, FILL_GHOST_CELLS);
}

// override virtual tagcells function for chi and phi gradients 
void MatterSFLevel::tagCells (IntVectSet& a_tags)
{
    CH_TIME("GRAMRLevel::tagCells");
    if ( m_verbosity ) pout () << "GRAMRLevel::tagCells " << m_level << endl;

    fillAllGhosts(); //We need filled ghost cells to calculate gradients etc

    // Create tags based on undivided gradient of phi and chi
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
<<<<<<< HEAD
            //At the moment only base on gradient chi and phi
            if ((mod_grad_fab(iv,c_chi)* m_dx >= m_p.regrid_threshold_chi)
            || (mod_grad_fab(iv,c_phi) * m_dx >= m_p.regrid_threshold_phi))
=======
            //At the moment only base on gradient chi/chi^2
            if ((mod_grad_fab(iv,c_chi) >= m_p.regrid_threshold_chi)
            || (mod_grad_fab(iv,c_phi) >= m_p.regrid_threshold_phi))
>>>>>>> e35718f7204c50f5d173a1f8ceafde6e6e79f65a
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
