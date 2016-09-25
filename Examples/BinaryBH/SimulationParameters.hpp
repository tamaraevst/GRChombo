#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

//General includes
#include "ParmParse.H"

//Problem specific includes:
#include "CCZ4.hpp"
#include "BoostedBH.hpp"

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

        //Fill in BinaryBHParameters
        bh1_params.mass = massA;
        bh1_params.center = centerA;
        bh1_params.momentum = momentumA;
        bh2_params.mass = massB;
        bh2_params.center = centerB;
        bh2_params.momentum = momentumB;

        //Fill in he ccz4Parameters
        ccz4Params.kappa1 = kappa1;
        ccz4Params.kappa2 = kappa2;
        ccz4Params.kappa3 = kappa3;
        ccz4Params.shift_gamma_coeff = F;
        ccz4Params.lapse_advec_coeff = LapseAdvectionCoeff;
        ccz4Params.shift_advec_coeff = ShiftAdvectionCoeff;
        ccz4Params.beta_driver = eta;
    }

    //SimulationParameters.inc declares all variables and defines auto_read_params(ParmParse& pp)
#include "SimulationParameters.inc"

    //Collection of parameters necessary for the CCZ4 RHS
    CCZ4::params_t ccz4Params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
