#ifndef SIMULATIONPARAMETERS_HPP_
#error  "This file must only be included through SimulationParameters.hpp"
#endif
//TODO: this annoying piece of code will be auto-generated from a list of values, names and default-values. This way we also get proper handling of missing parameters

void auto_read_params(ParmParse& pp)
{
    pp.get("verbosity", verbosity);
    //Grid setup
    pp.get("L", L);
    pp.get("regrid_threshold", regrid_threshold);
    pp.get("regrid_threshold_K", regrid_threshold_K);
    pp.get("regrid_threshold_phi", regrid_threshold_phi);
    pp.getarr("isPeriodic", isPeriodic,0, SpaceDim);
    pp.get("num_ghosts", num_ghosts);
    pp.get("tag_buffer_size", tag_buffer_size);
    pp.get("dt_multiplier", dt_multiplier);

    // Lapse evolution
    pp.get("LapseAdvectionCoeff",LapseAdvectionCoeff);

    // Shift Evolution
    pp.get("ShiftBCoeff",ShiftBCoeff);
    pp.get("ShiftAdvectionCoeff",ShiftAdvectionCoeff);
    pp.get("F",F);
    pp.get("eta",eta);
    pp.get("SpatialBetaDriverRadius",SpatialBetaDriverRadius);

    //CCZ4 parameters
    pp.get("formulation", formulation);
    pp.get("kappa1",kappa1);
    pp.get("kappa2",kappa2);
    pp.get("kappa3",kappa3);
    pp.get("covariantZ4",covariantZ4);

    //Dissipation
    pp.get("sigma",sigma);

    //Initial data
    pp.get("G_Newton", G_Newton);
    pp.get("amplitudeSF", amplitudeSF);
    pp.getarr("centerSF", centerSF, 0, SpaceDim);
    pp.get("widthSF", widthSF);

    //Relaxation params
    pp.get("relaxtime", relaxtime);
    pp.get("relaxspeed", relaxspeed);

    //Misc
    pp.get("nan_check", nan_check);
    ignore_checkpoint_name_mismatch = false;
    pp.query("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch);
}

////// General parameters
int verbosity;
Real L; //Physical sidelength of the grid
int num_ghosts;
int tag_buffer_size; //Amount the tagged region is grown by
std::vector<bool> isPeriodic;
int nan_check;
bool ignore_checkpoint_name_mismatch; //ignore mismatch of variable names between restart file and program
Real dt_multiplier;

///// Problem specific parameters
Real regrid_threshold_K, regrid_threshold_phi, regrid_threshold;
//Lapse evolution
Real LapseAdvectionCoeff;
//ShiftEvolution
Real ShiftAdvectionCoeff, F, eta, SpatialBetaDriverRadius;
int ShiftBCoeff;
//CCZ4 parameters
int formulation;
Real kappa1, kappa2, kappa3;
int covariantZ4;
Real sigma; //Dissipation
//Initial data
double G_Newton;
Real amplitudeSF, widthSF;
std::vector<Real> centerSF;
//Relaxation params
Real relaxtime, relaxspeed;
