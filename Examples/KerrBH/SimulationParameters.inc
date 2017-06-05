// Last edited K Clough 16.10.16

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
    pp.get("lapse_advec_coeff", lapse_advec_coeff);
    pp.get("lapse_coeff", lapse_coeff);
    pp.get("lapse_power", lapse_power);

    // Shift Evolution
    pp.get("shift_advec_coeff", shift_advec_coeff);
    pp.get("shift_Gamma_coeff",shift_Gamma_coeff);
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

    //Initial Kerr data
    pp.get("kerr_mass", kerr_mass);
    pp.get("kerr_spin", kerr_spin);
    for(int i=0; i < CH_SPACEDIM; i++)
    {
        pp.get("kerr_center", kerr_center[i], i);
    }

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
Real lapse_power, lapse_coeff, lapse_advec_coeff;
//ShiftEvolution
Real shift_advec_coeff, shift_Gamma_coeff, eta, SpatialBetaDriverRadius;
//CCZ4 parameters
int formulation;
double kappa1, kappa2, kappa3;
int covariantZ4;
Real sigma; //Dissipation
//Initial data for Kerr
double G_Newton;
Real kerr_mass, kerr_spin;
std::array<double, CH_SPACEDIM> kerr_center;
