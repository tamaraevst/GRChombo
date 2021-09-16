/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{   c_chibg,
    c_xMom,
    c_BHMom,
    c_gaussbonnet,
    c_chernsimons,
    c_phianalytic,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
 "chi background", "xMom", "BHMom", "Gauss Bonnet", "Chern Simons", "Phi Analytic"};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
