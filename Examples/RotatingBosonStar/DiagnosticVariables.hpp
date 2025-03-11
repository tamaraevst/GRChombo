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

    c_rho,   // stress tensor components
    c_s1,
    c_s2,
    c_s3,
    c_s11,
    c_s12,
    c_s13,
    c_s22,
    c_s23,
    c_s33,

    c_Qphi_density,

    c_Fphi_flux,

    c_Sphi_source,

    c_weight1,
    c_weight2,

    c_profile1,
    c_profile2,

    c_testHam,

    c_dt_mod_phi,
    c_gamma_tt,

    c_mode_dtA_01_re,
    c_mode_dtA_01_im,
    c_mode_dtA_02_re,
    c_mode_dtA_02_im,
    c_mode_dtA_03_re,
    c_mode_dtA_03_im,
    c_mode_dtA_04_re,
    c_mode_dtA_04_im,
    c_mode_dtA_05_re,
    c_mode_dtA_05_im,
    c_mode_dtA_06_re,
    c_mode_dtA_06_im,
    c_mode_dtA_07_re,
    c_mode_dtA_07_im,
    c_mode_dtA_08_re,
    c_mode_dtA_08_im,
    c_mode_dtA_09_re,
    c_mode_dtA_09_im,
    c_mode_dtA_10_re,
    c_mode_dtA_10_im,
    c_mode_dtA_11_re,
    c_mode_dtA_11_im,
    c_mode_dtA_12_re,
    c_mode_dtA_12_im,
    c_mode_dtA_13_re,
    c_mode_dtA_13_im,
    c_mode_dtA_14_re,
    c_mode_dtA_14_im,
    c_mode_dtA_15_re,
    c_mode_dtA_15_im,
    c_mode_dtA_16_re,
    c_mode_dtA_16_im,
    c_mode_dtA_17_re,
    c_mode_dtA_17_im,
    c_mode_dtA_18_re,
    c_mode_dtA_18_im,
    c_mode_dtA_19_re,
    c_mode_dtA_19_im,
    c_mode_dtA_20_re,
    c_mode_dtA_20_im,

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

    "rho", "s1", "s2", "s3", "s11", "s12", "s13", "s22", "s23", "s33",

    "Qphi_density", "Fphi_flux", "Sphi_source",
    
    "weight1", "weight2",
        
    "profile1", "profile2",
    
    "testHam", 

    "dt_A_sq", "abs_gamma_tt",

    "mode_dtA_01_re", "mode_dtA_01_im", "mode_dtA_02_re", "mode_dtA_02_im", "mode_dtA_03_re", "mode_dtA_03_im", "mode_dtA_04_re", "mode_dtA_04_im", "mode_dtA_05_re", "mode_dtA_05_im", "mode_dtA_06_re", "mode_dtA_06_im", "mode_dtA_07_re", "mode_dtA_07_im", "mode_dtA_08_re", "mode_dtA_08_im", "mode_dtA_09_re", "mode_dtA_09_im", "mode_dtA_10_re", "mode_dtA_10_im", "mode_dtA_11_re", "mode_dtA_11_im", "mode_dtA_12_re", "mode_dtA_12_im", "mode_dtA_13_re", "mode_dtA_13_im", "mode_dtA_14_re", "mode_dtA_14_im", "mode_dtA_15_re", "mode_dtA_15_im", "mode_dtA_16_re", "mode_dtA_16_im", "mode_dtA_17_re", "mode_dtA_17_im", "mode_dtA_18_re", "mode_dtA_18_im", "mode_dtA_19_re", "mode_dtA_19_im", "mode_dtA_20_re", "mode_dtA_20_im"
    
     };
    
}

#endif /* DIAGNOSTICVARIABLES_HPP */
