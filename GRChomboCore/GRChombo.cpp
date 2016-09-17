#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>

using std::endl;

#include "GRChombo.hpp"
#include "BRMeshRefine.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "ParmParse.H"
#include "CH_Timer.H"

#include "FArrayBox.H"
#include "LevelData.H"
#include "LayoutIterator.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FourthOrderFillPatch.H"
#include "LevelFluxRegister.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "LevelRK4.H"
#include "Copier.H"
#include "computeNorm.H"


#include "FABView.H"

#include "FABDriver.hpp"
#include "Constraints.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "EnforceTfA.hpp"
#include "ComputeModGrad.hpp"
#include "NanCheck.hpp"

#ifdef USE_PAPI
#include "PapiProfilingInfo.hpp"
#endif

#include <fenv.h>
#if defined(__i386__) && defined(__SSE__)
#include <xmmintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "UsingNamespace.H"

#include "BoostedBH.hpp"
#include "BinaryBH.hpp"

#warning: GRChombo.cpp and hpp will be completely refactored and split into classes GRAMRLevel and CCZ4AMRLevel

/// Global variables for handling output:
static const char* pgmname = "GRChombo" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static int verbose = 0 ;

const int
GRChombo::s_num_comps;

const char*
GRChombo::s_state_names[s_num_comps] =
{
    "chi",

    "h11",
    "h12",
    "h13",
    "h22",
    "h23",
    "h33",

    "K",

    "A11",
    "A12",
    "A13",
    "A22",
    "A23",
    "A33",

    "Theta",

    "Gamma1",
    "Gamma2",
    "Gamma3",

    "lapse",

    "shift1",
    "shift2",
    "shift3",

    "B1",
    "B2",
    "B3",

    "Ham",
    "Mom1",
    "Mom2",
    "Mom3"
};

const int
GRChombo::s_num_ghosts;

GRChombo::GRChombo (const SimParams &a_p, int a_tagBufferSize, ProfilingInfo * a_profilingInfo)
: m_tagBufferSize(a_tagBufferSize), m_p(a_p), m_profilingInfo(a_profilingInfo)
{
    if ( verbose ) pout () << "GRChombo default constructor" << endl;
}

GRChombo::~GRChombo ()
{
}

void
GRChombo::define (AMRLevel* a_coarser_level_ptr,
                  const Box& a_problem_domain,
                  int a_level,
                  int a_ref_ratio)
{
    ProblemDomain physdomain(a_problem_domain);

    define(a_coarser_level_ptr, physdomain, a_level, a_ref_ratio);
}


void
GRChombo::define (AMRLevel* a_coarser_level_ptr,
                  const ProblemDomain& a_problem_domain,
                  int a_level,
                  int a_ref_ratio)
{
    if ( verbose ) pout () << "GRChombo::define " << a_level << endl;

    AMRLevel::define (a_coarser_level_ptr,
                      a_problem_domain,
                      a_level,
                      a_ref_ratio);

    if (a_coarser_level_ptr) {
        GRChombo *coarser_level_ptr = dynamic_cast<GRChombo *>(a_coarser_level_ptr);
        if (coarser_level_ptr) {
            m_dx = coarser_level_ptr->m_dx / Real(coarser_level_ptr->m_ref_ratio);
        } else {
            MayDay::Error ("in GRChombo::define: a_coarser_level_ptr is not castable to GRChombo*");
        }
    } else {
        m_dx = m_p.L / (a_problem_domain.domainBox().longside ());
    }
}


// advance by one timestep
Real
GRChombo::advance ()
{
    CH_TIME("GRChombo::advance");

    //if ( verbose )
    pout () << "GRChombo::advance " << m_level << " at " << m_time << endl;

    m_state_new.copyTo ( m_state_new.interval (),
                         m_state_old,
                         m_state_old.interval () );

    // The level classes take flux-register parameters, use dummy ones here
    LevelFluxRegister dummy_fr;
    LevelFluxRegister* coarser_fr = &dummy_fr;
    LevelFluxRegister* finer_fr   = &dummy_fr;

    // Undefined leveldata in case we need it
    const LevelData<FArrayBox> dummy_data;

    // Set arguments to dummy values and then fix if real values are available
    const LevelData<FArrayBox>* coarser_data_old = &dummy_data;
    const LevelData<FArrayBox>* coarser_data_new = &dummy_data;

    Real t_coarser_old = 0.0;
    Real t_coarser_new = 0.0;

    // A coarser level exists
    if (m_coarser_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_coarser_level_ptr);
        if (grchombo_ptr != NULL)
        {
            coarser_data_old = &grchombo_ptr->m_state_old;
            coarser_data_new = &grchombo_ptr->m_state_new;

            t_coarser_new = grchombo_ptr->m_time;
            t_coarser_old = t_coarser_new - grchombo_ptr->m_dt;
        }
        else
        {
            MayDay::Error ("in GRChombo::advance: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }

    if (m_finer_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_finer_level_ptr);
        if (grchombo_ptr != NULL)
        {
            RK4LevelAdvance(m_state_new, m_state_old,
                            grchombo_ptr->m_patcher.getTimeInterpolator(),
                            *coarser_data_old, t_coarser_old,
                            *coarser_data_new, t_coarser_new,
                            *coarser_fr, *finer_fr,
                            m_time, m_dt, *this);
        }
        else
        {
            MayDay::Error ("in GRChombo::advance: m_finer_level_ptr is not castable to GRChombo*");
        }

    } else {
        RK4LevelAdvance(m_state_new, m_state_old,
                        *coarser_data_old, t_coarser_old,
                        *coarser_data_new, t_coarser_new,
                        *coarser_fr, *finer_fr,
                        m_time, m_dt, *this);
    }

    // Now we need to fix the algebraic constraints
    const DisjointBoxLayout& level_domain = m_state_new.disjointBoxLayout();
    DataIterator dit0 = level_domain.dataIterator();
    int nbox = dit0.size();

    //Print nBox to give information on load balancing
    pout () << "Number of level " << m_level << " boxes on this rank: " << nbox << "." << endl;

    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(m_state_new, m_state_new, true);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(m_state_new, m_state_new, true);

    //Check for nan's
    if (m_p.nan_check) FABDriver<NanCheck>().execute(m_state_new, m_state_new, true, disable_simd());

    m_time += m_dt;
    return m_dt;

}


// methods used with LevelRK4
// evaluate d(soln)/dt at current time based on soln
void
GRChombo::evalRHS(TSoln& rhs, // d(soln)/dt based on soln
                  TSoln& soln, // soln at current time
                  TFR& fineFR,  // flux register w/ finer level
                  TFR& crseFR,  // flux register w/ crse level
                  const TSoln& oldCrseSoln, // old-time crse solution
                  Real oldCrseTime,    // old crse time
                  const TSoln& newCrseSoln,  // new-time crse solution
                  Real newCrseTime,   // new crse time
                  Real time,   // current time centering of soln
                  Real fluxWeight // weight to apply to fluxRegister updates
                 )
{
    CH_TIME("GRChombo::evalRHS");

    soln.exchange(m_exchange_copier);

    if (oldCrseSoln.isDefined()) {
        // Fraction "a_time" falls between the old and the new coarse times
        Real alpha = (time - oldCrseTime) / (newCrseTime - oldCrseTime);

        // Truncate the fraction to the range [0,1] to remove floating-point
        // subtraction roundoff effects
        Real eps = 0.04 * m_dt / m_ref_ratio;

        if (Abs(alpha) < eps)
        {
            alpha = 0.0;
        }

        if (Abs(1.0-alpha) < eps)
        {
            alpha = 1.0;
        }

        // Current time before old coarse time
        if (alpha < 0.0)
        {
            MayDay::Error( "GRChombo::evalRHS: alpha < 0.0");
        }

        // Current time after new coarse time
        if (alpha > 1.0)
        {
            MayDay::Error( "GRChombo::evalRHS: alpha > 1.0");
        }

        // Interpolate ghost cells from next coarser level using both space
        // and time interpolation
        m_patcher.fillInterp(soln,
                             alpha,
                             0,0,s_num_comps);
    }
    //Time and count the RHS if possible
    if (m_profilingInfo != NULL) m_profilingInfo->resetCounters();

    //Enforce the trace free alpha condition //TODO: Make sure it's enough to enforce it here! Then delete other place
    //where it's enforced
    FABDriver<EnforceTfA>().execute(soln, soln, true);

    //Enforce positive chi and alpha
    FABDriver<PositiveChiAndAlpha>().execute(soln, soln, true);

    //Calculate CCZ4 right hand side
    FABDriver<CCZ4>(m_p.ccz4Params, m_dx, m_p.sigma).execute(soln, rhs, false);

    if (m_profilingInfo != NULL) m_profilingInfo->readCounters();
}


// implements soln += dt*rhs
void
GRChombo::updateODE(TSoln& soln,
                    const TSoln& rhs,
                    Real dt)
{
    CH_TIME("GRChombo::updateODE");

    DataIterator dit0 = soln.dataIterator();
    int nbox = dit0.size();

    for(int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        FArrayBox& soln_fab = soln[di];
        const FArrayBox& rhs_fab = rhs[di];
        soln_fab.plus(rhs_fab, dt);
    }

    //Enforce the trace free alpha condition
    FABDriver<EnforceTfA>().execute(soln, soln, true);
}


// define data holder newSoln based on existingSoln,
// including ghost cell specification
void
GRChombo::defineSolnData(TSoln& newSoln,
                         const TSoln& existingSoln)
{
    newSoln.define (existingSoln.disjointBoxLayout(), existingSoln.nComp(),
                    existingSoln.ghostVect());
}


// define data holder for RHS based on existingSoln
// including ghost cell specification
// (which in most cases is no ghost cells)
void
GRChombo::defineRHSData(TSoln& newRHS, const TSoln& existingSoln)
{
    newRHS.define (existingSoln.disjointBoxLayout(), existingSoln.nComp());
}


/// copy data in src into dest
void
GRChombo::copySolnData(TSoln& dest, const TSoln& src)
{
    src.copyTo ( src.interval (),
                 dest,
                 dest.interval () );
}


// things to do after a timestep
void
GRChombo::postTimeStep ()
{
    if ( verbose ) pout () << "GRChombo::postTimeStep " << m_level << endl;

    if (m_finer_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_finer_level_ptr);
        if ( grchombo_ptr != NULL)
        {
            grchombo_ptr->m_coarse_average.averageToCoarse (m_state_new,
                                                            grchombo_ptr->m_state_new);
        }
        else
        {
            MayDay::Error ("in GRChombo::postTimeStep: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }
    if ( verbose ) pout () << "GRChombo::postTimeStep " << m_level << " finished" << endl;
}


// create tags
void
GRChombo::tagCells (IntVectSet& a_tags)
{
    CH_TIME("GRChombo::tagCells");

    if ( verbose ) pout () << "GRChombo::tagCells " << m_level << endl;

    // Since tags are calculated using only current time step data, use
    // the same tagging function for initialization and for regridding.

    // Create tags based on undivided gradient of phi
    IntVectSet local_tags;

    const DisjointBoxLayout& level_domain = m_state_new.disjointBoxLayout();
    // If there is a coarser level then interpolate undefined ghost cells
    if (m_coarser_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_coarser_level_ptr);
        if (grchombo_ptr != NULL)
        {
            FourthOrderFillPatch fofp(level_domain,
                                      grchombo_ptr->m_state_new.disjointBoxLayout(),
                                      s_num_comps,
                                      grchombo_ptr->m_problem_domain,
                                      grchombo_ptr->m_ref_ratio,
                                      s_num_ghosts);

            fofp.fillInterp(m_state_new,
                            grchombo_ptr->m_state_new,
                            0,
                            0,
                            s_num_comps);
        }
        else
        {
            MayDay::Error ("in GRChombo::tagCells: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }

    m_state_new.exchange(Interval(0,s_num_comps-1));

    // Compute the constraint monitor
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
        for (int z = zmin; z <= zmax; ++z) {
            for (int y = ymin; y <= ymax; ++y) {
                for (int x = xmin; x <= xmax; ++x)
                {
                    IntVect iv(x,y,z);
                    //At the moment only base on gradient chi/chi^2
                    if (mod_grad_fab(iv,c_chi)/pow(state_fab(iv,c_chi),2) >= m_p.regridmax)
                    {
                        // local_tags |= is not thread safe.
#pragma omp critical
                        {
                            local_tags |= iv;
                        }
                    }
                }//x
            }//y
        }//z
    }

    local_tags.grow(m_tagBufferSize);

    // Need to do this in two steps unless a IntVectSet::operator &=
    // (ProblemDomain) operator is defined
    Box local_tags_box = local_tags.minBox();
    local_tags_box &= m_problem_domain;
    local_tags &= local_tags_box;

    a_tags = local_tags;
}

// create tags at initialization

void
GRChombo::tagCellsInit (IntVectSet& a_tags)
{
    tagCells(a_tags);
}

DisjointBoxLayout
GRChombo::loadBalance(const Vector<Box>& a_grids)
{
    CH_TIME("GRChombo::loadBalance");

    // load balance and create boxlayout
    Vector<int> procMap;

    // appears to be faster for all procs to do the loadbalance (ndk)
    LoadBalance(procMap,a_grids);

    if (verbose)
    {
        pout() << "GRChombo::::loadBalance: procesor map: " << endl;

        for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
            pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }

        pout() << endl;
    }

    DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
    dbl.close();

    return dbl;
}

// regrid

void
GRChombo::regrid (const Vector<Box>& a_new_grids)
{
    CH_TIME("GRChombo::regrid");

    if ( verbose ) pout () << "GRChombo::regrid " << m_level << endl;

    m_level_grids = a_new_grids;

    mortonOrdering(m_level_grids);
    const DisjointBoxLayout level_domain = m_grids = loadBalance (a_new_grids);

    for (LayoutIterator lit = level_domain.layoutIterator (); lit.ok (); ++lit)
    {
        if ( verbose ) pout () << level_domain[lit ()] << endl;
    }

    // save data for later copy
    m_state_new.copyTo(m_state_new.interval(),
                       m_state_old,
                       m_state_old.interval());

    // reshape state with new grids
    IntVect iv_ghosts = s_num_ghosts*IntVect::Unit;
    m_state_new.define (level_domain, s_num_comps, iv_ghosts);

    // maintain interlevel stuff
    m_exchange_copier.exchangeDefine(level_domain, iv_ghosts);
    m_coarse_average.define (level_domain,
                             s_num_comps,
                             m_ref_ratio);
    m_fine_interp.define (level_domain,
                          s_num_comps,
                          m_ref_ratio,
                          m_problem_domain);

    if (m_coarser_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_coarser_level_ptr);
        if (grchombo_ptr != NULL)
        {
            m_patcher.define(level_domain,
                             grchombo_ptr->m_grids,
                             s_num_comps,
                             grchombo_ptr->problemDomain(),
                             m_ref_ratio,
                             s_num_ghosts);

            // interpolate from coarser level
            m_fine_interp.interpToFine (m_state_new,
                                        grchombo_ptr->m_state_new);
        }
        else
        {
            MayDay::Error ("in GRChombo::regrid: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }
    // copy from old state
    m_state_old.copyTo (m_state_old.interval (),
                        m_state_new,
                        m_state_new.interval () );
    m_state_old.define (level_domain, s_num_comps, iv_ghosts);
}


// initialize grid
void
GRChombo::initialGrid (const Vector<Box>& a_new_grids)
{
    CH_TIME("GRChombo::initialGrid");

    if (verbose) pout () << "GRChombo::initialGrid " << m_level << endl;

    m_level_grids = a_new_grids;

    const DisjointBoxLayout level_domain = m_grids = loadBalance (a_new_grids);

    for (LayoutIterator lit = level_domain.layoutIterator (); lit.ok (); ++lit)
    {
        if ( verbose ) pout () << level_domain[lit ()] << endl;
    }

    IntVect iv_ghosts = s_num_ghosts*IntVect::Unit;
    m_state_new.define (level_domain, s_num_comps, iv_ghosts);
    m_state_old.define (level_domain, s_num_comps, iv_ghosts);

    m_exchange_copier.exchangeDefine(level_domain, iv_ghosts);
    m_coarse_average.define (level_domain,
                             s_num_comps,
                             m_ref_ratio);
    m_fine_interp.define (level_domain,
                          s_num_comps,
                          m_ref_ratio,
                          m_problem_domain);

    if (m_coarser_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_coarser_level_ptr);
        if (grchombo_ptr != NULL)
        {
            m_patcher.define(level_domain,
                             grchombo_ptr->m_grids,
                             s_num_comps,
                             grchombo_ptr->problemDomain(),
                             m_ref_ratio,
                             s_num_ghosts);
        }
        else
        {
            MayDay::Error ("in GRChombo::initialGrid: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }
}

// initialize data
// this is where the initial data for the metric, curvature etc are input
void
GRChombo::initialData ()
{
    CH_TIME("GRChombo::initialData");

    if (verbose) pout () << "GRChombo::initialData " << m_level << endl;

#warning: All the code below will go ... it should be moved to a compute class

    DataIterator dit0 = m_state_new.dataIterator();
    int nbox = dit0.size();
    //#pragma omp parallel for default(shared) schedule(guided)
    for(int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex di = dit0[ibox];
        FArrayBox& state_fab = m_state_new[di];
        Box b = state_fab.box ();
        const IntVect& smallEnd = b.smallEnd();
        const IntVect& bigEnd = b.bigEnd();

        const int xmin = smallEnd[0];
        const int ymin = smallEnd[1];
        const int zmin = smallEnd[2];

        const int xmax = bigEnd[0];
        const int ymax = bigEnd[1];
        const int zmax = bigEnd[2];

#pragma omp parallel for collapse(3) schedule(static) default(shared)
        for (int ivz = zmin; ivz <= zmax; ++ivz) {
            for (int ivy = ymin; ivy <= ymax; ++ivy) {
                for (int ivx = xmin; ivx <= xmax; ++ivx)
                {
                    IntVect iv(ivx,ivy,ivz);

                    for (int comp = 0; comp < m_state_new.nComp (); ++comp)
                    {
                        state_fab (iv,comp) = 0;
                    }

                    //Note: Cell centred!
                    Real x = (iv[0] + 0.5) * m_dx;
                    Real y = (iv[1] + 0.5) * m_dx;
                    Real z = (iv[2] + 0.5) * m_dx;

                    BoostedBH bh1(m_p.massA, m_p.centerA, m_p.momentumA);
                    BoostedBH bh2(m_p.massB, m_p.centerB, m_p.momentumB);

                    BinaryBH<BoostedBH> binary(bh1, bh2);

                    // Metric conformal factor
                    const Real psi = binary.psi(x, y, z);
                    state_fab(iv, c_chi) = pow(psi, -4);

                    // Conformal metric is flat
                    state_fab(iv,c_h11) = 1;
                    state_fab(iv,c_h22) = 1;
                    state_fab(iv,c_h33) = 1;

                    // Maximal slicing
                    state_fab(iv,c_K) = 0;

                    // Extrinsic curvature
                    Real Aij[3][3] = {{0}};
                    Real BOOK2BSSN = pow(psi, -6);
                    binary.Aij(x, y, z, Aij);


                    state_fab(iv,c_A11) = BOOK2BSSN * Aij[0][0];
                    state_fab(iv,c_A12) = BOOK2BSSN * Aij[0][1];
                    state_fab(iv,c_A13) = BOOK2BSSN * Aij[0][2];
                    state_fab(iv,c_A22) = BOOK2BSSN * Aij[1][1];
                    state_fab(iv,c_A23) = BOOK2BSSN * Aij[1][2];
                    state_fab(iv,c_A33) = BOOK2BSSN * Aij[2][2];

                    // Lapse
                    state_fab(iv,c_lapse) = 1;
                }//ivx
            }//ivy
        }//ivz
    }
}


// things to do after initialization
void
GRChombo::postInitialize ()
{
    if ( verbose ) pout () << "GRChombo::postInitialize " << m_level << endl;
}


// write checkpoint header
#ifdef CH_USE_HDF5
void
GRChombo::writeCheckpointHeader (HDF5Handle& a_handle) const
{
    if ( verbose ) pout () << "GRChombo::writeCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.m_int ["num_components"] = s_num_comps;
    char comp_str[30];
    for (int comp = 0; comp < s_num_comps; ++comp)
    {
        sprintf (comp_str, "component_%d", comp);
        header.m_string[comp_str] = s_state_names[comp];
    }
    header.writeToFile(a_handle);

    if ( verbose ) pout () << header << endl;
}

void
GRChombo::writeCheckpointLevel (HDF5Handle& a_handle) const
{
    CH_TIME("GRChombo::writeCheckpointLevel");

    if ( verbose ) pout () << "GRChombo::writeCheckpointLevel" << endl;

    char level_str[20];
    sprintf (level_str, "%d", m_level);
    const std::string label = std::string ("level_") + level_str;

    a_handle.setGroup (label);

    HDF5HeaderData header;

    header.m_int  ["ref_ratio"]   = m_ref_ratio;
    header.m_int  ["tag_buffer_size"] = m_tagBufferSize;
    header.m_real ["dx"]          = m_dx;
    header.m_real ["dt"]          = m_dt;
    header.m_real ["time"]        = m_time;
    header.m_box  ["prob_domain"] = m_problem_domain.domainBox();

    // Setup the periodicity info
    D_TERM6(if (m_problem_domain.isPeriodic(0))
            {
            header.m_int ["is_periodic_0"] = 1;
            }
            else
            {
            header.m_int ["is_periodic_0"] = 0;
            } ,

            if (m_problem_domain.isPeriodic(1))
            {
            header.m_int ["is_periodic_1"] = 1;
            }
            else
            {
            header.m_int ["is_periodic_1"] = 0;
            } ,

            if (m_problem_domain.isPeriodic(2))
            {
            header.m_int ["is_periodic_2"] = 1;
            }
            else
            {
                header.m_int ["is_periodic_2"] = 0;
            } ,

                if (m_problem_domain.isPeriodic(3))
                {
                    header.m_int ["is_periodic_3"] = 1;
                }
                else
                {
                    header.m_int ["is_periodic_3"] = 0;
                } ,

                    if (m_problem_domain.isPeriodic(4))
                    {
                        header.m_int ["is_periodic_4"] = 1;
                    }
                    else
                    {
                        header.m_int ["is_periodic_4"] = 0;
                    } ,

                        if (m_problem_domain.isPeriodic(5))
                        {
                            header.m_int ["is_periodic_5"] = 1;
                        }
                        else
                        {
                            header.m_int ["is_periodic_5"] = 0;
                        } );

    header.writeToFile(a_handle);

    if ( verbose ) pout () << header << endl;

    write (a_handle, m_state_new.boxLayout ());
    write (a_handle, m_state_new, "data");
}

void
GRChombo::preCheckpointLevel ()
{
#warning should fill ghost cells here (put into convenient call)
    FABDriver<Constraints>(m_dx).execute(m_state_new, m_state_new, false);
}


void
GRChombo::readCheckpointHeader  (HDF5Handle& a_handle)
{
    CH_TIME("GRChombo::readCheckpointHeader");

    if ( verbose ) pout () << "GRChombo::readCheckpointHeader" << endl;

    HDF5HeaderData header;
    header.readFromFile(a_handle);

    if ( verbose ) pout () << "hdf5 header data:" << endl;
    if ( verbose ) pout () << header << endl;

    // read number of components
    if (header.m_int.find ("num_components") == header.m_int.end())
    {
        MayDay::Error ("GRChombo::readCheckpointHeader: checkpoint file does not have num_components");
    }
    int num_comps = header.m_int ["num_components"];
    if (num_comps != s_num_comps)
    {
        MayDay::Error ("GRChombo::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

    // read component names
    std::string state_name;
    char comp_str[60];
    for (int comp = 0; comp < s_num_comps; ++comp)
    {
        sprintf (comp_str, "component_%d", comp);
        if (header.m_string.find (comp_str) == header.m_string.end())
        {
            MayDay::Error ("GRChombo::readCheckpointHeader: checkpoint file does not have enough component names");
        }
        state_name = header.m_string [comp_str];
        if (state_name != s_state_names[comp])
        {
//            if (m_p.ignoreNameMismatch)
//                MayDay::Warning
//            else
                MayDay::Error("GRChombo::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }

}

void
GRChombo::readCheckpointLevel (HDF5Handle& a_handle)
{
    CH_TIME("GRChombo::readCheckpointLevel");

    if ( verbose ) pout () << "GRChombo::readCheckpointLevel" << endl;

    //  char* level_str = new char[int (log10 (m_level+1))+1];
    char level_str[20];
    sprintf (level_str, "%d", m_level);
    const std::string label = std::string ("level_") + level_str;
    //  delete[] level_str;

    a_handle.setGroup (label);

    HDF5HeaderData header;
    header.readFromFile (a_handle);

    if ( verbose ) pout () << "hdf5 header data:" << endl;
    if ( verbose ) pout () << header << endl;

    // read refinement ratio
    if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain ref_ratio");
    }
    m_ref_ratio = header.m_int ["ref_ratio"];

    if ( verbose ) pout () << "read ref_ratio = " << m_ref_ratio << endl;

    // read dx
    if (header.m_real.find("dx") == header.m_real.end())
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain dx");
    }
    m_dx = header.m_real ["dx"];

    if ( verbose ) pout () << "read dx = " << m_dx << endl;

    // read dt
    //if (header.m_real.find("dt") == header.m_real.end())
    //{
    //  MayDay::Error("GRChombo::readCheckpointLevel: file does not contain dt");
    //}
    //m_dt = header.m_real ["dt"];
    //Since we have fixed time steping it is better to take this from parameter file:
    computeInitialDt();

    if ( verbose ) pout () << "read dt = " << m_dt << endl;

    // read time
    if (header.m_real.find("time") == header.m_real.end())
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain time");
    }
    m_time = header.m_real ["time"];

    if ( verbose ) pout () << "read time = " << m_time << endl;

    // read problem domain
    if (header.m_box.find("prob_domain") == header.m_box.end())
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain prob_domain");
    }
    Box domainBox = header.m_box ["prob_domain"];

    // Get the periodicity info -- this is more complicated than it really
    // needs to be in order to preserve backward compatibility
    bool isPeriodic[SpaceDim];
    D_TERM6(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
            {
            isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
            }
            else
            {
            isPeriodic[0] = false;
            } ,

            if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
            {
            isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
            }
            else
            {
            isPeriodic[1] = false;
            } ,

            if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
            {
            isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
            }
            else
            {
                isPeriodic[2] = false;
            } ,

                if (!(header.m_int.find("is_periodic_3") == header.m_int.end()))
                {
                    isPeriodic[3] =  (header.m_int["is_periodic_3"] == 1);
                }
                else
                {
                    isPeriodic[3] = false;
                } ,

                    if (!(header.m_int.find("is_periodic_4") == header.m_int.end()))
                    {
                        isPeriodic[4] =  (header.m_int["is_periodic_4"] == 1);
                    }
                    else
                    {
                        isPeriodic[4] = false;
                    } ,

                        if (!(header.m_int.find("is_periodic_5") == header.m_int.end()))
                        {
                            isPeriodic[5] =  (header.m_int["is_periodic_5"] == 1);
                        }
                        else
                        {
                            isPeriodic[5] = false;
                        } );

    m_problem_domain = ProblemDomain(domainBox, isPeriodic);

    // read grids
    Vector<Box> grids;
    const int grid_status = read (a_handle, grids);
    if (grid_status != 0)
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain a Vector<Box>");
    }

    // create level domain
    const DisjointBoxLayout level_domain = m_grids = loadBalance (grids);

    if ( verbose ) pout () << "read level domain: " << endl;
    LayoutIterator lit = level_domain.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
        const Box& b = level_domain[lit()];
        if ( verbose ) pout () << lit().intCode() << ": " << b << endl;
        m_level_grids.push_back(b);
    }
    if ( verbose ) pout () << endl;

    // maintain interlevel stuff
    IntVect iv_ghosts = s_num_ghosts*IntVect::Unit;
    m_exchange_copier.exchangeDefine(level_domain, iv_ghosts);
    m_coarse_average.define (level_domain,
                             s_num_comps,
                             m_ref_ratio);
    m_fine_interp.define (level_domain,
                          s_num_comps,
                          m_ref_ratio,
                          m_problem_domain);

    if (m_coarser_level_ptr != NULL)
    {
        GRChombo* grchombo_ptr =
            dynamic_cast<GRChombo*> (m_coarser_level_ptr);
        if (grchombo_ptr != NULL)
        {
            m_patcher.define(level_domain,
                             grchombo_ptr->m_grids,
                             s_num_comps,
                             grchombo_ptr->problemDomain(),
                             m_ref_ratio,
                             s_num_ghosts);
        }
        else
        {
            MayDay::Error ("in GRChombo::readCheckpointLevel: m_coarser_level_ptr is not castable to GRChombo*");
        }
    }

    // reshape state with new grids
    m_state_new.define (level_domain, s_num_comps, iv_ghosts);
    const int data_status = read<FArrayBox> (a_handle,
                                             m_state_new,
                                             "data",
                                             level_domain);
    if (data_status != 0)
    {
        MayDay::Error("GRChombo::readCheckpointLevel: file does not contain state data");
    }
    m_state_old.define (level_domain, s_num_comps, iv_ghosts);
}


void
GRChombo::writePlotLevel (HDF5Handle& a_handle) const
{
    if ( verbose ) pout () << "GRChombo::writePlotLevel" << endl;
}


void
GRChombo::writePlotHeader (HDF5Handle& a_handle) const
{
    if ( verbose ) pout () << "GRChombo::writePlotHeader" << endl;
}
#endif


// compute dt
Real
GRChombo::computeDt ()
{
    if ( verbose ) pout () << "GRChombo::computeDt " << m_level << endl;
    return m_dt;
}


// compute dt with initial data
Real
GRChombo::computeInitialDt ()
{
    if ( verbose ) pout () << "GRChombo::computeInitialDt " << m_level << endl;

    m_dt = m_initial_dt_multiplier * m_dx;
    return m_dt;
}


class GRChomboFactory : public AMRLevelFactory
{
public:
    GRChomboFactory(ProfilingInfo * profilingInfo=NULL);

    virtual AMRLevel* new_amrlevel() const;

    virtual
        ~GRChomboFactory();

protected:
    Real m_dt_multiplier;
    SimParams m_p;
    int m_tagBufferSize;
    ProfilingInfo* m_profilingInfo;
};


GRChomboFactory::GRChomboFactory (ProfilingInfo * a_profilingInfo):
    m_profilingInfo (a_profilingInfo)
{
    ParmParse pp;
    pp.get("dt_multiplier", m_dt_multiplier);
    m_p.readParams(pp);
    pp.get("tag_buffer_size", m_tagBufferSize);
}
GRChomboFactory::~GRChomboFactory ()
{
}


// "virtual constructor"
AMRLevel*
GRChomboFactory::new_amrlevel() const
{
    GRChombo*
        grchombo_ptr = new GRChombo (m_p, m_tagBufferSize, m_profilingInfo);
    grchombo_ptr->initialDtMultiplier(m_dt_multiplier);
    return (static_cast <AMRLevel*> (grchombo_ptr));
}


/// Prototypes:
int
runGRChombo();

int
main(int argc ,char* argv[])
{
#if defined(FE_NOMASK_ENV) && !defined(__INTEL_COMPILER)
    fesetenv(FE_NOMASK_ENV);
    fedisableexcept(/* FE_OVERFLOW | */ FE_UNDERFLOW | FE_INEXACT);
#elif defined(__i386__) && defined(__SSE__)
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~(_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO));
#endif

#ifdef CH_MPI
    // Start MPI
    MPI_Init(&argc,&argv);
#ifdef CH_AIX
    H5dont_atexit();
#endif
    // setChomboMPIErrorHandler();
#endif

    int rank, number_procs;

#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

    if (rank == 0)
    {
        pout() << " number_procs = " << number_procs << endl;
#ifdef _OPENMP
        pout() << " threads = " << omp_get_max_threads() << endl;
#endif
    }

    if (argc < 2)
    {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    pp.get("verbose", verbose);
    if ( verbose )
        pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

    int status = runGRChombo();

    if ( status == 0 )
        pout() << indent << pgmname << " finished." << endl ;
    else
        pout() << indent << pgmname << " failed with return code " << status << endl ;

#ifdef CH_MPI
    // Exit MPI
    dumpmemoryatexit();
    MPI_Finalize();
#endif

    return status ;
}


int
runGRChombo ()
{
    ParmParse pp;

    IntVect ivN = IntVect::Unit;
    D_TERM6(
            int N1; pp.get("N1", N1);
            ivN[0] = N1-1;,
            int N2; pp.get("N2", N2);
            ivN[1] = N2-1;,
            int N3; pp.get("N3", N3);
            ivN[2] = N3-1;,
            int N4; pp.get("N4", N4);
            ivN[3] = N4-1;,
            int N5; pp.get("N5", N5);
            ivN[4] = N5-1;,
            int N6; pp.get("N6", N6);
            ivN[5] = N6-1;)

        Box problem_domain (IntVect::Zero,
                            ivN);

    ProblemDomain physdomain(problem_domain);

    // set periodic in all directions
    std::vector<bool> isPeriodic;
    pp.getarr("isPeriodic", isPeriodic, 0, SpaceDim);
    for (int dir=0; dir<SpaceDim; dir++)
    {
        physdomain.setPeriodic(dir, isPeriodic[dir]);
    }

    int max_level;
    pp.get("max_level", max_level);

    Vector<int> ref_ratios;
    pp.getarr("ref_ratio",ref_ratios,0,max_level+1);

#ifdef USE_PAPI
    int numEvents = 5;
    int events[5] = {PAPI_L1_DCM,PAPI_L1_ICM,PAPI_L2_DCM,PAPI_L2_ICM,PAPI_L3_TCM};
    PapiProfilingInfo papiProfilingInfo(events,numEvents);
    papiProfilingInfo.startCounters();
    GRChomboFactory grchombo_fact(&papiProfilingInfo);
#else
    GRChomboFactory grchombo_fact;
#endif


    AMR amr;
    amr.define(max_level, ref_ratios, physdomain, &grchombo_fact);

    // To preserve proper nesting we need to know the maximum ref_ratio.
    int max_ref_ratio = ref_ratios[0];
    for (int i = 1; i < max_level+1; ++i) {
        max_ref_ratio = std::max(max_ref_ratio, ref_ratios[i]);
    }

    int grid_buffer_size = std::ceil(
                                     ((double) GRChombo::s_num_ghosts)/ (double) max_ref_ratio
                                    ) + 3 /* stencil dist. */;
    amr.gridBufferSize(grid_buffer_size);

    int checkpoint_interval;
    pp.get("checkpoint_interval", checkpoint_interval);
    amr.checkpointInterval(checkpoint_interval);

    // Number of coarse time steps from one regridding to the next
    Vector<int> regrid_intervals;
    pp.getarr("regrid_interval",regrid_intervals,0,max_level+1);
    amr.regridIntervals(regrid_intervals);

    if (pp.contains("max_grid_size"))
    {
        int max_grid_size;
        pp.query("max_grid_size", max_grid_size);
        amr.maxGridSize(max_grid_size);
    }

    if (pp.contains("block_factor"))
    {
        int block_factor;
        pp.query("block_factor", block_factor);
        amr.blockFactor(block_factor);
    }

    if (pp.contains("fill_ratio"))
    {
        Real fill_ratio;
        pp.query("fill_ratio", fill_ratio);
        amr.fillRatio(fill_ratio);
    }

    if (pp.contains("max_dt_grow"))
    {
        Real max_dt_grow;
        pp.query("max_dt_grow", max_dt_grow);
        amr.maxDtGrow(max_dt_grow);
    }

    if (pp.contains("dt_tolerance_factor"))
    {
        Real dt_tolerance_factor;
        pp.query("dt_tolerance_factor", dt_tolerance_factor);
        amr.dtToleranceFactor(dt_tolerance_factor);
    }

    if (pp.contains("chk_prefix"))
    {
        std::string prefix;
        pp.query("chk_prefix",prefix);
        amr.checkpointPrefix(prefix);
    }

    amr.verbosity(verbose);

    // Set up input files
    if (!pp.contains("restart_file")) {
        amr.setupForNewAMRRun();
    } else {
        std::string restart_file;
        pp.query("restart_file",restart_file);

#ifdef CH_USE_HDF5
        HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        amr.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("GRChombo restart only defined with hdf5");
#endif
    }

    Real stop_time;
    pp.get("stop_time", stop_time);

    int max_steps;
    pp.get("max_steps", max_steps);

    amr.run(stop_time, max_steps);

    amr.conclude();

#ifdef USE_PAPI
    //Read one final time
    papiProfilingInfo.readShutdown();
#endif

    return 0 ;
}

