/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This file calculates CCZ4 geometric quantities (or a similar 3+1 split).
#ifndef CCZ4GEOMETRYMODIFIEDGR_HPP_
#define CCZ4GEOMETRYMODIFIEDGR_HPP_

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"
#include "Tensor.hpp"
#include "CCZ4Geometry.hpp"

/*This class has functions to assist with the computation of Chern Simons 
and Gauss Bonnet scalars*/
class CCZ4GeometryModifiedGR
{   
    public:

    //struct for some derivative expressions of the lapse
    template <class data_t> struct lapse_t
    {   
        Tensor<2, data_t> covd2lapse;
        data_t tr_covd2lapse;      
        Tensor<2, data_t> tr_free_covd2lapse;             
    };

    //struct for \delta_t expressions in the decomposition of the CS and GB terms
    template <class data_t> struct evolution_t
    {
        Tensor<2, data_t> A;
        data_t K;
        data_t chi;
    };

    /* This fuction computes D_k A_ij in conformal variables. Note that the corresponding quantity for it is covd_Aij[i][j][k].
    Note the order of brackets for the indices [i][j][k] corrresponding to the prder of _k _ij in the actual tensor!
    */
    template <class data_t, template <typename> class vars_t>
    static Tensor<3, data_t>
    compute_covd_A(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        Tensor<3, data_t> covd_A;

        using namespace TensorAlgebra;

        FOR3(i, j, k)
        {
            //This menans D_k A_{ij}
            covd_A[i][j][k] = (1.0 / (2.0 * vars.chi * vars.chi)) * (2.0 * vars.chi * d1.A[i][j][k] + 
                vars.A[j][k] * d1.chi[i] + vars.A[i][k] * d1.chi[j]);

            FOR2(m, n)
            {
                covd_A[i][j][k] += (1.0 / vars.chi) * ( - chris.ULL[m][i][k] * vars.A[m][j] - chris.ULL[m][j][k] * vars.A[i][m]) -
                            (h_UU[n][m] / (2.0 * vars.chi * vars.chi)) * (vars.A[m][j] * vars.h[i][k] * d1.chi[n] + 
                            vars.A[i][m] * vars.h[j][k] * d1.chi[n]);
            }    
        }
        return covd_A;
    }
    /* This function computes some of the useful quantities related to the lapse, more precisely we have:
    1) covd2lapse = D_i D_j \alpha 
    2) tr_covd2lapse = D_i D^i \alpha
    3) tr_free_covd2lapse = [D_i D_j \alpha]^{TF}, i.e. the trace free part of 1)
    */
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static lapse_t<data_t> 
    compute_covd2lapse_quantities(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {   
        using namespace TensorAlgebra;

        lapse_t<data_t> out;

        data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

        Tensor<2, data_t> covdtilde2lapse;

        //First, expression for covd2lapse
        FOR2(k, l)
        {
            covdtilde2lapse[k][l] = d2.lapse[k][l];
            FOR1(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
            out.covd2lapse[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
        }

        //Next, expression for tr_covd2lapse
        out.tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
        FOR1(i)
        {
            out.tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
            FOR1(j)
            {
                out.tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
            }
        }

        //Finally, expression for tr_free_covd2lapse
        FOR2(i, j)
        {
            out.tr_free_covd2lapse[i][j] = out.covd2lapse[i][j];
        }
        make_trace_free(out.tr_free_covd2lapse, vars.h, h_UU);
        return out;
    }
    
    /* There are some evolution terms in Gauss-Bonnet and Chern-Simons scalars, i.e. \partial_t A_ij, \partial_t K and \partial_t K
    We present their BSSN decomposition in the function below. Note that advection terms are substracted from them in the formula.
    Hence, we include this substraction into our definiton of "evolution terms". More precisely, we have teh following variables calculated in this function:
    1) A = \partial_t A_{ij} - \beta^s \delta_s A_{ij}
    2) K = \partial_t K - \beta^s \delta_s K
    3) chi = \partial_t \chi - \beta^s \delta_s \chi
    */
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static evolution_t<data_t>
    rhs_evolution_quantities(const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
                    const diff2_vars_t<Tensor<2, data_t>> &d2, const chris_t<data_t> &chris)
    {
    using namespace TensorAlgebra;

    evolution_t<data_t> rhs;

    //Calling some functions needed for computation
    auto h_UU = compute_inverse_sym(vars.h);

    Tensor<1, data_t> Z0 = 0.0;
    auto ricci =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);

    auto lapse_derivatives = compute_covd2lapse_quantities(vars, d1, d2, h_UU, chris); 

    data_t divshift = compute_trace(d1.shift);

    Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    // A^{ij} A_{ij} = A^2 quantity 
    data_t tr_A2 = compute_trace(vars.A, A_UU);

    //First, term chi.
    rhs.chi = (2.0 / GR_SPACEDIM) * vars.chi * (vars.lapse * vars.K - divshift);

    //Next, term A. It has a trace free term, which we compute below and call Adot_TF
    Tensor<2, data_t> Adot_TF;
    FOR2(i, j)
    {
        Adot_TF[i][j] =
            -lapse_derivatives.covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, vars.h, h_UU);

    FOR2(i, j)
    {
        rhs.A[i][j] = Adot_TF[i][j] +
                      vars.A[i][j] * (vars.lapse * vars.K -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR1(k)
        {
            rhs.A[i][j] +=
                vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i];
            FOR1(l)
            {
                rhs.A[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
            }
        }
    }

    //Finally, term K.
    rhs.K = vars.lapse * (tr_A2 + vars.K * vars.K /  GR_SPACEDIM) -
                lapse_derivatives.tr_covd2lapse;

    return rhs;
}

    //This function compoutes electric term, E_ij, present in the Chern Simons scalar
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    Tensor<2, data_t>
    compute_chern_simons_electric_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU)
    {
        using namespace TensorAlgebra;

        Tensor<2, data_t> K_tensor;

        auto chris = compute_christoffel(d1.h, h_UU);

        //Calling some functions needed for computation
        // auto covd2lapse = compute_covd2lapse_quantities(vars, d1, d2, h_UU);
        
        // auto rhs = rhs_evolution_quantities(vars, d1, d2);

        Tensor<2, data_t> E = 0.0; //the Chern-Simons electric term

        Tensor<1, data_t> Z0 = 0.0;
        auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);
        
        // Tensor<2, data_t> ricciscalar_TF;
        // FOR2(i, j)
        //  {
        //     ricciscalar_TF[i][j] = ricci.LL[i][j];
        //  }
        // make_trace_free(ricciscalar_TF, vars.h, h_UU);

        Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);
        data_t tr_A2 = compute_trace(vars.A, A_UU); // A^{ij} A_{ij}

        //Finally, compute the electric part of the Chern-Simons term E_{ij}.
        // FOR2(i, j)
        // {
        //     Eij[i][j] = (1.0 / (2.0 * vars.chi * vars.chi * vars.lapse)) * (- rhs.chi * vars.A[i][j] + vars.chi * rhs.A[i][j]) + 
        //         (1.0 / 2.0) * (covd2lapse.tr_free_covd2lapse[i][j] * (1.0 / vars.lapse) + ricciscalar_TF[i][j]) +
        //         (1.0 / (3.0 * vars.chi)) * vars.h[i][j] * tr_A2 + (1.0 / (6 * vars.chi)) * vars.K * vars.A[i][j];
        //     FOR1(k)
        //     {
        //         Eij[i][j] += (1.0 / (2.0 * vars.chi * vars.lapse)) * (- vars.A[k][j] * d1.shift[k][i] - vars.A[i][k] * d1.shift[k][j]);
        //     }
        // }

        FOR2(i, j)
        {
            K_tensor[i][j] = vars.A[i][j] / vars.chi +
                         1. / 3. * (vars.h[i][j] * vars.K) / vars.chi;
        }
        FOR2(i, j)
        {
            E[i][j] = ricci.LL[i][j] + vars.K * K_tensor[i][j];
            
            FOR2(k, l)
            {
                E[i][j] +=
                     -K_tensor[i][k] * K_tensor[l][j] * h_UU[k][l] * vars.chi;
            }
        }

        TensorAlgebra::make_trace_free(E, vars.h, h_UU);

        return E;
    }
    
    //Fucntion for computing the magnetic term, B_ij, in the Chern-Simons scalar. The same B_ij is then encountered in Gauss-Bonnet. 
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    Tensor<2, data_t>
    compute_magnetic_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU)

    {
        using namespace TensorAlgebra;

        Tensor<2, data_t> K_tensor;
        Tensor<3, data_t> d1_K_tensor;
        Tensor<3, data_t> covariant_deriv_K_tensor;

        Tensor<2, data_t> B; 

        auto chris = compute_christoffel(d1.h, h_UU);
        const Tensor<3, data_t> chris_phys =
        compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

        const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);
        
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

        FOR4(i, j, k, l)
        {
            B[i][j] +=
                epsilon3_LUU[i][k][l] * covariant_deriv_K_tensor[l][j][k];
        }

        TensorAlgebra::make_symmetric(B);

        // const auto covd_A = compute_covd_A(vars, d1, h_UU);

        //Compute the magnetic term here
        // FOR2(i, j)
        // {
        //     FOR2(k, l)
        //     {
        //         Bij[i][j] = (1.0 / 2.0) * (epsilon3_LUU[i][k][l] * covd_A[j][l][k] + epsilon3_LUU[j][k][l] * covd_A[i][l][k]);
        //     }  
        // }
        return B;
    }
    
    /* This function is for computing almost all of the decomposed Gauss Bonnet (GB) term, exclusing M and B terms.
    It is introduced to do most of the tedious computations seperately.
    */
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    data_t GB_scalar(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU)
{
     
        using namespace TensorAlgebra;

        data_t GB_scalar;

        auto chris = compute_christoffel(d1.h, h_UU);

        auto lapse_quantities = compute_covd2lapse_quantities(vars, d1, d2, h_UU, chris);

        Tensor<1, data_t> Z0 = 0.0;
        auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);
        
        Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);
        data_t tr_A2 = compute_trace(vars.A, A_UU);

        auto rhs = rhs_evolution_quantities(vars, d1, d2, chris);

        //This part is the first line of eq (2.4) in the note, the term that is multipled with H^{GR}.
        GB_scalar = 4.0 / 3.0 * (ricci.scalar - tr_A2 + 2.0 / 3.0 * vars.K * vars.K) *
            (1 / vars.lapse * (rhs.K  + lapse_quantities.tr_covd2lapse) - tr_A2 - 1.0 / 3.0 * vars.K * vars.K);

        //This is the second line of eq (2.4) in the note, the term that is multiplied with E^{GR}.
        FOR4(i, j, k, l)
        {
            GB_scalar += 8 * (vars.chi * vars.chi * h_UU[i][k] * h_UU[j][l] * ricci.LL[i][j] - 
            vars.chi * A_UU[i][k] * h_UU[l][j] * vars.A[j][i] + 1.0 / 3.0 * vars.chi * (vars.K * A_UU[k][l] +
            h_UU[k][l] * tr_A2)) * ( - vars.A[k][l] / (vars.chi * vars.chi) * rhs.chi + 
            1.0 / vars.chi * rhs.A[k][l]  + 1.0 / vars.chi * (vars.A[i][l] * d1.shift[i][k] - 
            vars.A[i][k] * d1.shift[i][l]) + 1.0 / (vars.lapse) * lapse_quantities.tr_free_covd2lapse[k][l] + 
            1.0 / vars.chi * vars.A[j][k] * h_UU[j][i] * vars.A[i][l]);
        }

        //Finally, we add B and M terms that are on the third line of eq (2.4) of the note.
        Tensor<3, data_t> covd_A = compute_covd_A(vars, d1, h_UU, chris);

        //First, define B
        const auto B_ij = compute_magnetic_term(vars, d1, d2, h_UU);
        
        //Raise index for B^ij and convert to conformal variables
        Tensor<2, data_t> B_UU;

        FOR4(i, j, k ,l)
        {
            B_UU[i][j] = vars.chi * vars.chi * h_UU[k][i] * h_UU[l][j] * B_ij[k][l];
        }
        
        //Next define M (consisting of GR momentum constraint)
        Tensor<1, data_t> M;

        //This is M_k
        FOR3(i, k, s)
        {
            M[k] = vars.chi * h_UU[i][s] * covd_A[k][i][s] - 2.0 / 3.0 * d1.K[k];
        }
        //Raise index M^k
        Tensor<1, data_t> M_U;
        FOR2(k, s)
        {
            M_U[k] = vars.chi * h_UU[k][s] * M[k];
        }
        
        //Finally, add M and B contributions!
        FOR4(i, j, k, l)
        {
            GB_scalar += - 8 * B_ij[i][j] * B_UU[i][j] + 4 * M[i] * M_U[i];
        }

        return GB_scalar;
    }

///This is copied from Weyl4 code; the function is protected there for calculating \epsilon_i^jk
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

#endif /* CCZ4GEOMETRYMODIFIEDGR_HPP_ */