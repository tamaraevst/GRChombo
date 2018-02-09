/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4VARS_HPP_
#define CCZ4VARS_HPP_

#include "ADMVars.hpp"
#include "BSSNVars.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

namespace CCZ4Vars
{
/// CCZ4 variables
/** This struct collects all the CCZ4 variables. It's main use is to make a
 *local, nicely laid-out, copy of the CCZ4 variables for the current grid
 *cell (Otherwise, this data would only exist on the grid in the huge,
 *flattened Chombo array). To this end, CCZ4::Vars inherits from VarsBase
 *which contains functionality to connect the local copy of the variables
 *with values in the Chombo grid.
 **/
template <class data_t>
struct VarsNoGauge : public BSSNVars::VarsNoGauge<data_t>
{
    // Additional variable is Theta
    data_t Theta; //!< CCZ4 quantity associated to hamiltonian constraint

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        BSSNVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_Theta, Theta);
    }
};

template <class data_t>
struct VarsWithGauge : public BSSNVars::VarsWithGauge<data_t>
{
    // Additional variable is Theta
    data_t Theta; //!< CCZ4 quantity associated to hamiltonian constraint

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        BSSNVars::VarsWithGauge<data_t>::enum_mapping(mapping_function);
        define_enum_mapping(mapping_function, c_Theta, Theta);
    }
};

/// 2nd derivatives are only calculated for a small subset defined by
/// Deriv2Vars
template <class data_t>
struct Diff2VarsNoGauge : public ADMVars::Diff2VarsNoGauge<data_t>
{
};

template <class data_t>
struct Diff2VarsWithGauge : public ADMVars::Diff2VarsWithGauge<data_t>
{
};
}

#endif /* CCZ4VARS_HPP_ */
