//This file incldues several functions that need to be called to
//set up the runs but aren't very interesting for the normal user.

#include "parstream.H" //Gives us pout()
#include <iostream>
using std::endl;
using std::cerr;
#include "ParmParse.H"
#include "AMR.H"
#include "AMRLevelFactory.H"


//TODO (MK): There is a lot of clutter still in this file ... get rid of everything that's not necessary
//and comment on the rest (so that we don't carry around useless code forever

void mainSetup(int argc ,char* argv[])
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
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);
}

void mainFinalize()
{
#ifdef CH_MPI
    // Exit MPI
    dumpmemoryatexit();
    MPI_Finalize();
#endif
}

void setupAMRObject(AMR& amr, AMRLevelFactory& a_factory)
{
    //Some hard-coded parameters:
    //The buffer is width of ghost cells + additional_grid_buffer
    //and defines the minimum number of level l cells there have to be
    //between level l+1 and level l-1
    const int additional_grid_buffer = 3;

    ParmParse pp;

    IntVect ivN = IntVect::Unit;
    // Setup the grid size
    for (int dir=0; dir<SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf (dir_str, "N%d", dir);
        int N;
        pp.get(dir_str, N);
        ivN[dir] = N-1;
    }

    Box problem_domain (IntVect::Zero, ivN);
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

    //Define the AMR object
    amr.define(max_level, ref_ratios, physdomain, &a_factory);

    // To preserve proper nesting we need to know the maximum ref_ratio.
    int max_ref_ratio = ref_ratios[0];
    for (int i = 1; i < max_level+1; ++i)
    {
        max_ref_ratio = std::max(max_ref_ratio, ref_ratios[i]);
    }

    int grid_buffer_size = std::ceil(
                                     ((double) GRChombo::s_num_ghosts)/ (double) max_ref_ratio
                                    ) + additional_grid_buffer;
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

    int verbosity;
    pp.get("verbosity", verbosity);
    amr.verbosity(verbosity);

    // Set up input files
    if (!pp.contains("restart_file"))
    {
        amr.setupForNewAMRRun();
    }
    else
    {
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
}
