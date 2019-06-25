/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "BinaryBH.hpp"
#include "BoxLoops.hpp"
#include "CCZ4.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "Constraints.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SetValue.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"
#include "ConstraintViolations.hpp"

void BinaryBHLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBHLevel::initialData " << m_level << endl;

    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::preCheckpointLevel()
{
    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side and set constraints to zero to avoid
    // undefined values
    BoxLoops::loop(
        make_compute_pack(CCZ4(m_p.ccz4_params, m_dx, m_p.sigma),
                          SetValue(0, Interval(c_Ham, NUM_VARS - 1))),
        a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                            const FArrayBox &current_state)
{
    BoxLoops::loop(
        ChiExtractionTaggingCriterion(m_dx, m_level, m_p.extraction_params),
        current_state, tagging_criterion);
}

void BinaryBHLevel::doAnalysis()
{
    CH_TIME("BinaryBHLevel::specificPostTimeStep");

    // Populate the Weyl Scalar values and Constraint violations on the grid
    fillAllGhosts();
    BoxLoops::loop(make_compute_pack(Constraints(m_dx),
                   Weyl4(m_p.extraction_params.extraction_center, m_dx)),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

    bool called_in_do_analysis = true;

    // Do the extraction on the min extraction level
    if (m_level == m_p.extraction_params.min_extraction_level
            && m_p.activate_extraction == 1)
    {

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        WeylExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                     m_restart_time, called_in_do_analysis);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    if (m_level == 0)
    {
        ConstraintViolations constraint_violations(c_Ham,
            Interval(c_Mom1, c_Mom3), &m_gr_amr, m_p.coarsest_dx, m_dt, m_time,
            m_restart_time, "ConstraintViolations.dat", called_in_do_analysis);
        constraint_violations.execute();
        auto violations = constraint_violations.get_norms();
        pout() << "L2 norms of constraint violations:\n";
        pout() << "Ham: " << violations.first << "\t";
        pout() << "Mom: " << violations.second << "\n";
    }

}

// Things to do before a plot level - need to calculate the Weyl scalars
void BinaryBHLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(Weyl4(m_p.extraction_params.extraction_center, m_dx),
                       m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    }
}

// Specify if you want any plot files to be written, with which vars
void BinaryBHLevel::specificWritePlotHeader(std::vector<int> &plot_states) const
{
    plot_states = {c_chi, c_Weyl4_Re, c_Weyl4_Im};
}
