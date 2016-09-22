#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

//General includes
#include "ParmParse.H"

//Problem specific includes:
#include "CCZ4.hpp"
#include "BoostedBH.hpp"

struct SimulationParameters
{
  void readParams(ParmParse &pp)
  {
#warning: this annoying piece of code will be auto-generated from a list of values, names and default-values. This way we also get proper handling of missing parameters
     //Grid setup
     pp.get("L", L);
     pp.get("regridmax", regridmax);
     pp.getarr("isPeriodic", isPeriodic,0, SpaceDim);
     pp.get("num_ghosts", num_ghosts);

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

     //Fill in BinaryBHParameters
     bh1_params.mass = massA;
     bh1_params.center = centerA;
     bh1_params.momentum = momentumA;
     bh2_params.mass = massB;
     bh2_params.center = centerB;
     bh2_params.momentum = momentumB;

     //Misc
     pp.get("nan_check", nan_check);
     pp.get("ignore_name_mismatch", ignore_name_mismatch);

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
  int num_ghosts;
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
  int nan_check;
  bool ignore_name_mismatch;

  //Collection of parameters necessary for the CCZ4 RHS
  CCZ4::params_t ccz4Params;
  BoostedBH::params_t bh2_params;
  BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
