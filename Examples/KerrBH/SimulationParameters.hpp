// Last edited K Clough 16.10.16

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

//General includes
#include "ParmParse.H"

//Problem specific includes:
#include "CCZ4.hpp"
#include "KerrBH.hpp"

class SimulationParameters
{
public:
    SimulationParameters(ParmParse& pp)
    {
        readParams(pp);
    }

    void readParams(ParmParse &pp)
    {
        //The automatically generated read parameters code defined in SimulationParameters.inc
        auto_read_params(pp);

        //Fill in KerrBH Parameters
        kerr_params.mass = kerr_mass;
        kerr_params.center = kerr_center;
        kerr_params.spin = kerr_spin;

        //Fill in the ccz4Parameters
        ccz4Params.kappa1 = kappa1;
        ccz4Params.kappa2 = kappa2;
        ccz4Params.kappa3 = kappa3;
        ccz4Params.shift_Gamma_coeff = shift_Gamma_coeff;
        ccz4Params.shift_advec_coeff = shift_advec_coeff;
        ccz4Params.eta = eta;
        ccz4Params.lapse_power = lapse_power;
        ccz4Params.lapse_coeff = lapse_coeff;
        ccz4Params.lapse_advec_coeff = lapse_advec_coeff;
    }

    //SimulationParameters.inc declares all variables and defines auto_read_params(ParmParse& pp)
#include "SimulationParameters.inc"

    //Collection of parameters necessary for the CCZ4 RHS
    CCZ4::params_t ccz4Params;
    KerrBH::params_t kerr_params;

};

#endif /* SIMULATIONPARAMETERS_HPP_ */
