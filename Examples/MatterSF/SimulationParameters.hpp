#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

//General includes
#include "ParmParse.H"

//Problem specific includes:
#include "CCZ4.hpp"
#include "BubbleSF.hpp"

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

        //Fill in the SFMatterParamters
        sfm_params.amplitudeSF = amplitudeSF;
        sfm_params.centerSF = centerSF;
        sfm_params.widthSF = widthSF;

        //Fill in the ccz4Parameters
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
    BubbleSF::params_t sfm_params;
  
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
