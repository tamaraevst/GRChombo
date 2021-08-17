/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This file calculates CCZ4 geometric quantities (or a similar 3+1 split).
#ifndef CCZ4GEOMETRYMODIFIEDGR_HPP_
#define CCZ4GEOMETRYMODIFIEDGR_HPP_

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "BSSNVars.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
// #include "CCZ4.hpp"

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

template <class data_t> struct modifiedscalar_t
{
    data_t starR_R;        // Pontryagin scalar
    data_t RGB;             // Gauss-Bonnet term
};

template <class data_t> struct covd2lapse_t
{   
    Tensor<2, data_t> covd2lapse;
    data_t tr_covd2lapse;      
    Tensor<2, data_t> tr_free_covd2lapse;             
};

template <class data_t> struct evolution_t
{
    Tensor<2, data_t> A;
    data_t K;
    data_t chi;
};

class CCZ4GeometryModifiedGR 
{  
    public:

    template <class data_t> using Vars = BSSNVars::VarsWithGauge<data_t>;

    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsWithGauge<data_t>;

    // template <class data_t> void compute(Cell<data_t> current_cell) const;
    
    template <class data_t, template <typename> class vars_t>
    static Tensor<3, data_t>
    compute_covd_Aij(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);

    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static covd2lapse_t<data_t> 
    compute_covd2lapse_quantities(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);
    
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static evolution_t<data_t>
    rhs_evolution_quantities(const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
                    const diff2_vars_t<Tensor<2, data_t>> &d2);

    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    Tensor<2, data_t>
    compute_chern_simons_electric_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);
    
    //For computing the magnetic part of the Chern-Simons term B_{ij}.
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    Tensor<2, data_t>
    compute_magnetic_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);
    

    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    data_t GB_scalar(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);

    // Function to compute modified scalars: Gauss Bonnet and Chern Simons
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static modifiedscalar_t<data_t>
    compute_modified_scalars(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris);

    // template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;
    template <class data_t, template <typename> class vars_t>
    Tensor<3, data_t> compute_epsilon3_LUU(const vars_t<data_t> &vars,
                                           const Tensor<2, data_t> &h_UU) const;
};

#include "CCZ4GeometryModifiedGR.impl.hpp"

#endif /* CCZ4GEOMETRYMODIFIEDGR_HPP_ */
