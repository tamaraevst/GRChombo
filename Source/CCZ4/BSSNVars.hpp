/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BSSNVARS_HPP_
#define BSSNVARS_HPP_

#include "ADMVars.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

namespace BSSNVars
{
/// BSSN variables
/** This struct collects all the BSSN variables. It's main use is to make a
 *local, nicely laid-out, copy of the BSSN variables for the current grid
 *cell (Otherwise, this data would only exist on the grid in the huge,
 *flattened Chombo array). To this end, BSSN::Vars inherits from VarsBase
 *which contains functionality to connect the local copy of the variables
 *with values in the Chombo grid.
 **/
template <class data_t> struct VarsNoGauge : public ADMVars::VarsNoGauge<data_t>
{
    Tensor<1, data_t> Gamma; //!< Conformal connection functions

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        ADMVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, GRInterval<c_Gamma1, c_Gamma3>(),
                            Gamma);
    }
};

template <class data_t> struct VarsWithGauge : public VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<1, data_t> B; //!< \f$B^i = \partial_t \beta^i\f$, this is used
                         //! for second order shift conditions

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        VarsNoGauge<data_t>::enum_mapping(mapping_function);
        // Scalars
        define_enum_mapping(mapping_function, c_lapse, lapse);

        // Vectors
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(),
                            shift);
        define_enum_mapping(mapping_function, GRInterval<c_B1, c_B3>(), B);
    }
};

/// 2nd derivatives are only calculated for a small subset defined by
/// Deriv2Vars
/** Making this split speeds up the code significantly */
template <class data_t>
struct Diff2VarsNoGauge : public ADMVars::Diff2VarsNoGauge<data_t>
{
};

template <class data_t>
struct Diff2VarsWithGauge : public ADMVars::Diff2VarsWithGauge<data_t>
{
};
}

#endif /* BSSNVARS_HPP_ */
