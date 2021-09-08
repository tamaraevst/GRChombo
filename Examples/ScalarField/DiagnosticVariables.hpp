/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom,

    c_chernsimons,

    c_gaussbonnet_1,

    c_gaussbonnet_2,

    c_phianalytic,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom",
    
    "ChernSimons",
    
    "GaussBonnet1",

    "GaussBonnet2",

    "Phi Analytic"

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
