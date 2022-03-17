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

    c_Mom1,
    c_Mom2,
    c_Mom3,

//    c_Moms_abs_terms,

//    c_Ham_abs_terms,
 
    c_chernsimons,

    c_gaussbonnet,

    c_phianalytic,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1", "Mom2", "Mom3",

//    "MomAbs",

//    "HamAbs",

    "chernsimons",
    
    "gaussbonnet",

    "phianalytic"

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
