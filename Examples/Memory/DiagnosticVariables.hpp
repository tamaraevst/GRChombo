/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef DIAGNOSTICVARIABLES_HPP
 #define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_mod_phi, // the complex scalar field modulus

    c_Madm,
    c_Jadm,

    c_N, // Noether Charge integrand

    c_Ham,

    c_Mom1,
    c_Mom2,
    c_Mom3,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_profile1,
    c_profile2,

    c_gxx,
    c_gxy,
    c_gxz,
    c_gyy, 
    c_gyz,
    c_gzz,

    c_gtxx,
    c_gtxy,
    c_gtxz,
    c_gtyy, 
    c_gtyz,
    c_gtzz,

    c_grxx,
    c_grxy,
    c_grxz,
    c_gryy, 
    c_gryz,
    c_grzz,

    c_shifttx,
    c_shiftty,
    c_shifttz,

    c_shiftrx,
    c_shiftry,
    c_shiftrz,

    c_lapset, 
    c_lapser,


    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {

    "mod_phi",

    "Madm",   "Jadm",

    "N",

    "Ham",    "Mom1",   "Mom2",   "Mom3",

    "Weyl4_Re",  "Weyl4_Im",

    "profile1", "profile2",

    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",

    "gtxx", "gtxy", "gtxz", "gtyy", "gtyz", "gtzz",

    "grxx", "grxy", "grxz", "gryy", "gryz", "grzz",

    "shifttx", "shiftty", "shifttz", 

    "shiftrx", "shiftry", "shiftrz",

    "lapset", "lapser",
    
     };
    
    }

#endif /* DIAGNOSTICVARIABLES_HPP */
