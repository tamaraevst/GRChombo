#ifndef SIMULATIONPARAMETERS_HPP_
#error  "This file must only be included through SimulationParameters.hpp"
#endif

void auto_read_params(GRParmParse& pp)
{
    pp.load("verbosity", verbosity);
    //Grid setup
    pp.load("L", L);
    pp.load("regrid_threshold", regrid_threshold);
    pp.load("regrid_threshold_K", regrid_threshold_K);
    pp.load("regrid_threshold_phi", regrid_threshold_phi);
    pp.load("isPeriodic", isPeriodic);
    pp.load("num_ghosts", num_ghosts);
    pp.load("tag_buffer_size", tag_buffer_size, 3);
    pp.load("dt_multiplier", dt_multiplier);

    // Lapse evolution
    pp.load("lapse_advec_coeff", lapse_advec_coeff);
    pp.load("lapse_coeff", lapse_coeff);
    pp.load("lapse_power", lapse_power);

    // Shift Evolution
    pp.load("shift_advec_coeff", shift_advec_coeff);
    pp.load("shift_Gamma_coeff",shift_Gamma_coeff);
    pp.load("eta", eta);
    pp.load("SpatialBetaDriverRadius",SpatialBetaDriverRadius);

    //CCZ4 parameters
    pp.load("formulation", formulation);
    pp.load("kappa1",kappa1);
    pp.load("kappa2",kappa2);
    pp.load("kappa3",kappa3);
    pp.load("covariantZ4",covariantZ4);

    //Dissipation
    pp.load("sigma",sigma);

    //Boson Star initial data params
    pp.load("central_amplitude_CSF", central_amplitude_CSF);
    pp.load("abs_error", abs_error);
    pp.load("rel_error", rel_error);
    pp.load("initial_step_size", initial_step_size);
    pp.load("max_radius", max_radius);
    pp.load("binary_search_tol", binary_search_tol);
    pp.load("max_binary_search_iter", max_binary_search_iter, 1e3);
    pp.load("star_centre", star_centre, {0.5*L, 0.5*L, 0.5*L});

    //Potential params
    pp.load("scalar_mass", scalar_mass);
    pp.load("phi4_coeff", phi4_coeff, 0.0);

    //Misc
    pp.load("nan_check", nan_check, 1);
    pp.load("ignore_checkpoint_name_mismatch", ignore_checkpoint_name_mismatch, false);
}

////// General parameters
int verbosity;
Real L; //Physical sidelength of the grid
int num_ghosts; // must be at least 3 for KO dissipation
int tag_buffer_size; //Amount the tagged region is grown by
std::array<bool, CH_SPACEDIM> isPeriodic;
int nan_check;
bool ignore_checkpoint_name_mismatch; //ignore mismatch of variable names between restart file and program
Real dt_multiplier;

///// Problem specific parameters
Real regrid_threshold_K, regrid_threshold_phi, regrid_threshold;
//Lapse evolution
Real lapse_power, lapse_coeff, lapse_advec_coeff;
//ShiftEvolution
Real shift_advec_coeff, shift_Gamma_coeff, eta, SpatialBetaDriverRadius;
//CCZ4 parameters
int formulation;
double kappa1, kappa2, kappa3;
int covariantZ4;
Real sigma; //Dissipation
//Initial data for matter and potential
double G_Newton;
double central_amplitude_CSF, abs_error, rel_error, initial_step_size;
double max_radius, binary_search_tol, max_binary_search_iter;
double scalar_mass, phi4_coeff;
std::array<double, CH_SPACEDIM> star_centre;
//Relaxation params
Real relaxtime, relaxspeed;
