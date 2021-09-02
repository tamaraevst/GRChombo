/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDSCALARSGEOMETRY_HPP
#define MODIFIEDSCALARSGEOMETRY_HPP

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"
#include "Tensor.hpp"
#include "CCZ4Geometry.hpp"

/*This class has functions to assist with the computation of Chern Simons 
and Gauss Bonnet scalars*/

/*The computation of E and B here is the same as in Weyl4 class implementation.
    This comes down to E and B being simplified to the same expressions as in the GR case. 
    Computation of this is shown in my note. Hence, we are indeed safe to use the same 
    expressions as in the Weyl4 class. While it is possible to inherit private functions of Weyl4 class, 
    I think it will make calling the class ComputeModifiedScalars using functions of this class
    too cubersome in terms of arguments. Furthermore, we will be using these scalars in ScalatField,
    so I would like to keep these functions public. Perhaps I am missing some better structure of the code.
    In practice, if we choose to have non-GR evolution equations and constraints, this
    class will have to be changed along with definitons of E and B. This is also why, I would
    like to leave these fucntions (which are geometrically significant to CS and GB decomposition) public*/

class ModifiedScalarsGeometry
{   
    public:

    // Struct for the gravito electric and magnetic parts E and B
    template <class data_t> struct EBterms_t   
    {
        Tensor<2, data_t> E; // gravito-electric part
        Tensor<2, data_t> B; // gravito-magnetic part
    };

    template <class data_t> struct modified_scalars_t
    {
        data_t GB_scalar;
        data_t CS_scalar;
    };

    public:
    //This function computes terms E and B present in the scalars decompositions
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    EBterms_t<data_t> compute_EB_terms(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2, const Tensor<3, data_t> &epsilon3_LUU,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        using namespace TensorAlgebra;
        EBterms_t<data_t> out;
        
        //Extrinsic curvature and derivatives
        Tensor<2, data_t> K_tensor;
        Tensor<3, data_t> d1_K_tensor;
        Tensor<3, data_t> covariant_deriv_K_tensor;

        const Tensor<3, data_t> chris_phys =
        compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);
        Tensor<1, data_t> Z0 = 0.0;
        auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);

        FOR(i, j)
        {
            out.E[i][j] = 0.0;
            out.B[i][j] = 0.0;
        }

        Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);
        data_t tr_A2 = compute_trace(vars.A, A_UU); // A^{ij} A_{ij}

        FOR2(i, j)
        {
            K_tensor[i][j] = vars.A[i][j] / vars.chi +
                         1. / 3. * (vars.h[i][j] * vars.K) / vars.chi;
            FOR(k)
            {
                d1_K_tensor[i][j][k] = d1.A[i][j][k] / vars.chi -
                                   d1.chi[k] / vars.chi * K_tensor[i][j] +
                                   1. / 3. * d1.h[i][j][k] * vars.K / vars.chi +
                                   1. / 3. * vars.h[i][j] * d1.K[k] / vars.chi;
            }
        }

        // covariant derivative of K
        FOR3(i, j, k)
        {
            covariant_deriv_K_tensor[i][j][k] = d1_K_tensor[i][j][k];

            FOR(l)
            {
                covariant_deriv_K_tensor[i][j][k] +=
                -chris_phys[l][k][i] * K_tensor[l][j] -
                chris_phys[l][k][j] * K_tensor[i][l];
            }
        }

        FOR2(i, j)
        {
            out.E[i][j] = ricci.LL[i][j] + vars.K * K_tensor[i][j];
            
            FOR2(k, l)
            {
                out.E[i][j] +=
                     -K_tensor[i][k] * K_tensor[l][j] * h_UU[k][l] * vars.chi;

                out.B[i][j] +=
                epsilon3_LUU[i][k][l] * covariant_deriv_K_tensor[l][j][k];
            }
        }

        make_trace_free(out.E, vars.h, h_UU);

        make_symmetric(out.B);

        return out;
    
    }
    
    /* This function is for computing GB scalar */
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    modified_scalars_t<data_t> mod_scalars(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
     
        using namespace TensorAlgebra;

        modified_scalars_t<data_t> out;

        out.GB_scalar = 0.0;
        out.CS_scalar = 0.0;

        const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

        auto ebterms = compute_EB_terms(vars, d1, d2, epsilon3_LUU, h_UU, chris);
        Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);
        data_t tr_A2 = compute_trace(vars.A, A_UU);

        FOR4(k, l, m , n)
        {
            out.GB_scalar += 8.0 * vars.chi * vars.chi * h_UU[k][m] * h_UU[l][n] * ebterms.E[m][n] * 
                        (ebterms.E[k][l] - (1.0)/(3.0 * vars.chi) * vars.h[k][l] * tr_A2) - 
                        8.0 * (vars.chi * vars.chi * h_UU[k][m] * h_UU[l][n] * ebterms.B[m][n] * ebterms.B[k][l]);

            out.CS_scalar = - 16.0 * vars.chi * vars.chi * h_UU[m][k] * h_UU[n][l] * ebterms.B[k][l] * ebterms.E[m][n];
        }
        return out;
    }

//This is 3D Levi Civita tensor (NB: not symbol), copied from Weyl4 code; the function is again protected
template <class data_t, template <typename> class vars_t>
Tensor<3, data_t> 
compute_epsilon3_LUU(const vars_t<data_t> &vars,
                    const Tensor<2, data_t> &h_UU) const
    {
        // raised normal vector, NB index 3 is time
        data_t n_U[4];
        n_U[3] = 1. / vars.lapse;
        FOR1(i) { n_U[i] = -vars.shift[i] / vars.lapse; }

        // 4D levi civita symbol and 3D levi civita tensor in LLL and LUU form
        const auto epsilon4 = TensorAlgebra::epsilon4D();
        Tensor<3, data_t> epsilon3_LLL;
        Tensor<3, data_t> epsilon3_LUU;

        // Projection of antisymmentric Tensor onto hypersurface - see 8.3.17,
        // Alcubierre
        FOR3(i, j, k)
        {
            epsilon3_LLL[i][j][k] = 0.0;
            epsilon3_LUU[i][j][k] = 0.0;
        }
        // projection of 4-antisymetric tensor to 3-tensor on hypersurface
        // note last index contracted as per footnote 86 pg 290 Alcubierre
        FOR3(i, j, k)
        {
            for (int l = 0; l < 4; ++l)
            {
                epsilon3_LLL[i][j][k] += n_U[l] * epsilon4[i][j][k][l] *
                                     vars.lapse / (vars.chi * sqrt(vars.chi));
            }
        }
    // rasing indices
        FOR3(i, j, k)
        {
            FOR2(m, n)
            {
                epsilon3_LUU[i][j][k] += epsilon3_LLL[i][m][n] * h_UU[m][j] *
                                        vars.chi * h_UU[n][k] * vars.chi;
            }
        }

    return epsilon3_LUU;
    }

};

#endif /* MODIFIEDSCALARSGEOMETRY_HPP_ */