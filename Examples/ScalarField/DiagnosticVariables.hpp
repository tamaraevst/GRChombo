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

    c_gaussbonnet,

    c_phianalytic,

    c_radius,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom",
    
    "ChernSimons",
    
    "GaussBonnet",

    "Phi Analytic",

    "Radius"

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
