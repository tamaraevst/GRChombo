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
    c_Ham_abs_terms,
    c_Mom,
    c_Moms_abs_terms,
 
    c_chernsimons,

    c_gaussbonnet,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "HamAbs",

    "Mom",

    "MomAbs",

    "chernsimons",
    
    "gaussbonnet"

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
