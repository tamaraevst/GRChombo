/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ADMVARS_HPP_
#define ADMVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

namespace ADMVars
{
/// ADM variables
/** This struct collects all the ADM variables. It's main use is to make a
 *local, nicely laid-out, copy of the ADM variables for the current grid
 *cell (Otherwise, this data would only exist on the grid in the huge,
 *flattened Chombo array). To this end, ADM::Vars inherits from VarsBase
 *which contains functionality to connect the local copy of the variables
 *with values in the Chombo grid.
 **/
template <class data_t> struct VarsNoGauge
{
    data_t chi;          //!< Conformal factor
    Tensor<2, data_t> h; //!< Conformal metric
    data_t K;            //!< Trace of the extrinsic curvature
    Tensor<2, data_t> A; //!< trace-free part of the rescale extrinsic
                         //! curvature, i.e. \f$\chi
                         //!(K_{ij})^{\mathrm{TF}}\f$

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        // Scalars
        define_enum_mapping(mapping_function, c_chi, chi);
        define_enum_mapping(mapping_function, c_K, K);

        // Symmetric 2-tensors
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h33>(), h);
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_A11, c_A33>(), A);
    }
};

template <class data_t> struct VarsWithGauge : public VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;

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
    }
};

/// 2nd derivatives are only calculated for a small subset defined by
/// Deriv2Vars
/** Making this split speeds up the code significantly */
template <class data_t> struct Diff2VarsNoGauge
{
    data_t chi;          //!< Conformal factor
    Tensor<2, data_t> h; //!< Conformal metric

    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        define_enum_mapping(mapping_function, c_chi, chi);
        define_symmetric_enum_mapping(mapping_function,
                                      GRInterval<c_h11, c_h33>(), h);
    }
};

template <class data_t>
struct Diff2VarsWithGauge : public Diff2VarsNoGauge<data_t>
{
    data_t lapse;
    Tensor<1, data_t> shift;

    /// Defines the mapping between members of Vars and Chombo grid
    /// variables (enum in User_Variables)
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        // Define the mapping from components of chombo grid to elements in
        // Vars. This allows to read/write data from the chombo grid into local
        // variables in Vars (which only exist for the current cell).

        using namespace VarsTools; // define_enum_mapping is part of VarsTools
        Diff2VarsNoGauge<data_t>::enum_mapping(mapping_function);
        // Scalars
        define_enum_mapping(mapping_function, c_lapse, lapse);

        // Vectors
        define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(),
                            shift);
    }
};
}

#endif /* ADMVARS_HPP_ */
