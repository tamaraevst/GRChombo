/* GRChombo Class
   */
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
#include "FourthOrderFineInterp.H"
#include "FourthOrderFillPatch.H"
#include "LevelFluxRegister.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "LevelRK4.H"
#include "Copier.H"
#include "computeNorm.H"

#include "ProfilingInfo.hpp"
#include "CCZ4.hpp"


#define NUMCOMP 25



struct SimParams
{
  void readParams(ParmParse &pp)
  {
     //Grid setup
     pp.get("L", L);
     pp.get("regridmax", regridmax);
     pp.getarr("isPeriodic", isPeriodic,0, SpaceDim);

     // Lapse evolution
     pp.get("LapseAdvectionCoeff",LapseAdvectionCoeff);

     // Shift Evolution
     pp.get("ShiftBCoeff",ShiftBCoeff);
     pp.get("ShiftAdvectionCoeff",ShiftAdvectionCoeff);
     pp.get("F",F);
     pp.get("eta",eta);
     pp.get("SpatialBetaDriverRadius",SpatialBetaDriverRadius);

     //CCZ4 parameters
     pp.get("kappa1",kappa1);
     pp.get("kappa2",kappa2);
     pp.get("kappa3",kappa3);
     pp.get("covariantZ4",covariantZ4);

     //Dissipation
     pp.get("sigma",sigma);

     //Initial data
     pp.get("massA", massA);
     pp.getarr("centerA", centerA, 0, SpaceDim);
     pp.getarr("momentumA", momentumA, 0, SpaceDim);
     pp.get("massB", massB);
     pp.getarr("centerB", centerB, 0, SpaceDim);
     pp.getarr("momentumB", momentumB, 0, SpaceDim);

     //Misc
     pp.get("nanCheck", nanCheck);

     //Fill in he ccz4Parameters
     ccz4Params.kappa1 = kappa1;
     ccz4Params.kappa2 = kappa2;
     ccz4Params.kappa3 = kappa3;
     ccz4Params.shift_gamma_coeff = F;
     ccz4Params.lapse_advec_coeff = LapseAdvectionCoeff;
     ccz4Params.shift_advec_coeff = ShiftAdvectionCoeff;
     ccz4Params.beta_driver = eta;
  }

  // Grid setup
  Real L, regridmax;
  std::vector<bool> isPeriodic;
  //Lapse evolution
  Real LapseAdvectionCoeff;
  //ShiftEvolution
  Real ShiftAdvectionCoeff, F, eta, SpatialBetaDriverRadius;
  int ShiftBCoeff;
  //CCZ4 parameters
  Real kappa1, kappa2, kappa3;
  int covariantZ4;
  //Dissipation
  Real sigma;
  //Initial data
  Real massA, massB;
  std::vector<Real> centerA, centerB;
  std::vector<Real> momentumA, momentumB;
  //Misc
  int nanCheck;

  //Collection of parameters necessary for the CCZ4 RHS
  CCZ4::params_t ccz4Params;
};

class GRChombo : public AMRLevel
{
public:
  GRChombo(const SimParams &a_p, int a_tagBufferSize, ProfilingInfo * a_profilingInfo = NULL);

  virtual
  ~GRChombo();

// define
  virtual
  void
  define(AMRLevel* a_coarser_level_ptr,
         const Box& a_problem_domain,
         int a_level,
         int a_ref_ratio);

// define
  virtual
  void
  define(AMRLevel* a_coarser_level_ptr,
         const ProblemDomain& a_problem_domain,
         int a_level,
         int a_ref_ratio);

// advance by one timestep
  virtual
  Real
  advance();

// things to do after a timestep
  virtual
  void
  postTimeStep();

// create tags
  virtual
  void
  tagCells(IntVectSet& a_tags) ;

// create tags at initialization
  virtual
  void
  tagCellsInit(IntVectSet& a_tags) ;

// regrid
  virtual
  void
  regrid(const Vector<Box>& a_new_grids);

// initialize grids
  virtual
  void
  initialGrid(const Vector<Box>& a_new_grids);

// initialize data
  virtual
  void
  initialData();

// things to do after initialization
  virtual
  void
  postInitialize();

#ifdef CH_USE_HDF5
  virtual
  void
  writeCheckpointHeader(HDF5Handle& a_handle) const;

  virtual
  void
  writeCheckpointLevel(HDF5Handle& a_handle) const;

  virtual
  void
  preCheckpointLevel();

  virtual
  void
  readCheckpointHeader(HDF5Handle& a_handle);

  virtual
  void
  readCheckpointLevel(HDF5Handle& a_handle);

  virtual
  void
  writePlotHeader(HDF5Handle& a_handle) const;

  virtual
  void
  writePlotLevel(HDF5Handle& a_handle) const;
#endif

// compute dt
  virtual
  Real
  computeDt();

// compute dt with initial data
  virtual
  Real
  computeInitialDt();

public:
// methods used with LevelRK4:
  typedef LevelData<FArrayBox> TSoln;
  typedef LevelFluxRegister TFR;

  // evaluate d(soln)/dt at current time based on soln
  void evalRHS(TSoln& rhs, // d(soln)/dt based on soln
               TSoln& soln, // soln at current time
               TFR& fineFR,  // flux register w/ finer level
               TFR& crseFR,  // flux register w/ crse level
               const TSoln& oldCrseSoln, // old-time crse solution
               Real oldCrseTime,    // old crse time
               const TSoln& newCrseSoln,  // new-time crse solution
               Real newCrseTime,   // new crse time
               Real time,   // current time centering of soln
               Real fluxWeight // weight to apply to fluxRegister updates
              );

  // implements soln += dt*rhs
  void updateODE(TSoln& soln,
                 const TSoln& rhs,
                 Real dt);

  // define data holder newSoln based on existingSoln,
  // including ghost cell specification
  void defineSolnData(TSoln& newSoln,
                      const TSoln& existingSoln);

  // define data holder for RHS based on existingSoln
  // including ghost cell specification
  // (which in most cases is no ghost cells)
  void defineRHSData(TSoln& newRHS, const TSoln& existingSoln);

  /// copy data in src into dest
  void copySolnData(TSoln& dest, const TSoln& src);

protected:
  DisjointBoxLayout
  loadBalance(const Vector<Box>& a_grids);

public:
// number of components of m_state
  static const int s_num_comps = NUMCOMP;
// names of components
  static const char* s_state_names[s_num_comps];
// number of ghost cells
  static const int s_num_ghosts = 3;
// periodicity information
  static const int s_periodicity[SpaceDim];

protected:
// state vector at old time
  LevelData<FArrayBox> m_state_old;
// state vector at new time
  LevelData<FArrayBox> m_state_new;
// grid spacing
  Real m_dx;
// tag buffer size
  int m_tagBufferSize;

// params
  SimParams m_p;

//Profiling info
  ProfilingInfo * m_profilingInfo;

// exchange copier
  Copier m_exchange_copier;

// interpolator for filling in ghost cells from the next coarser level
  FourthOrderFillPatch m_patcher;

  CoarseAverage m_coarse_average;
  FourthOrderFineInterp m_fine_interp;
  DisjointBoxLayout m_grids;
};
