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
#include "Weyl4.hpp"

template <class data_t> struct modfiedscalar_t
{
    data_t starR_R;        // Pontryagin scalar
    data_t RGB;             // Gauss-Bonnet term
};

class CCZ4GeometryModifiedGR
{   
    template <class data_t, template <typename> class vars_t>
    static Tensor<3, data_t>
    compute_covd_Aij(const vars_t<data_t> &vars,
                    const vars_t<Tensor<1, data_t>> &d1,
                    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        Tensor<3, data_t> covd_Aij = 0.0;
        FOR3(i, j, k)
        {
            //This menans D_k A_{ij}
            covd_Aij[i][j][k] = (1.0 / (2.0 * vars.chi * vars.chi)) * (2.0 * vars.chi * d1.A[i][j][k] + 
                + vars.A[k][j] * d1.chi[i] + vars.A[k][i] * d1.chi[j]);

            FOR1(m)
            {
                covd_Aij[i][j][k] -= (1.0 / vars.chi) * (chris.ULL[m][i][k] * vars.A[m][j] +
                    chris.ULL[m][j][k] * vars.A[i][m]);
                
                FOR1(n)
                {
                    covd_Aij[i][j][k] -= (h_UU[m][n] / (2.0 * vars.chi * vars.chi)) * (vars.A[n][m] * vars.h[i][k] * d1.chi[m] + 
                        vars.A[n][i] * vars.h[j][k] * d1.chi[m]);
                }
            }
        }
        return covd_Aij;
    }

    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static Tensor<2, data_t>
    compute_chern_simons_electric_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        using namespace TensorAlgebra;

        Tensor<2, data_t> Eij = 0.0; //the Chern-Simons electric term

        Tensor<1, data_t> Z0 = 0.;
        auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);

        data_t divshift = compute_trace(d1.shift);
        data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

        Tensor<2, data_t> covdtilde2lapse;
        Tensor<2, data_t> covd2lapse;
        FOR2(k, l)
        {
            covdtilde2lapse[k][l] = d2.lapse[k][l];
            FOR1(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
            covd2lapse[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
        }

        data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
        FOR1(i)
        {
            tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
            FOR1(j)
            {
                tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
            }
        }

        // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
        Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);
         data_t tr_A2 = compute_trace(vars.A, A_UU);

        //Compute the time evolution of chi. Note that its advection term is taken out, see formula for E_{ij}.
        data_t chi = (2.0 / GR_SPACEDIM) * vars.chi * (vars.lapse * vars.K - divshift);

        //Compute the time evolution of \tilde{A}_{ij}. Note that its advection term is taken out, see formula for E_{ij}.
        //Below the computation is split into traceless and trace parts.
        Tensor<2, data_t> Adot_TF;
        FOR2(i, j)
        {
            Adot_TF[i][j] =
                -covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
        }
        make_trace_free(Adot_TF, vars.h, h_UU);

        Tensor<2, data_t> TildeAij;
        FOR2(i, j)
        {
            TildeAij[i][j] = Adot_TF[i][j] +
                      vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) -
                                      (2.0 / GR_SPACEDIM) * divshift);
            FOR1(k)
            {
                TildeAij[i][j] +=
                vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i];
                FOR1(l)
                {
                    TildeAij[i][j] -=
                    2 * vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
                }
            }
        }

        //Finally, compute the electric part of the Chern-Simons term E_{ij}.
        FOR2(i, j)
        {
            Eij[i][j] = (1.0 / (2.0 * vars.chi * vars.chi)) * (chi * vars.A[i][j] - vars.chi * TildeAij) - 
                (2.0 * vars.chi * vars.chi) * ((tr_covd2lapse / vars.lapse) + ricci) +
                (1.0 / (3.0 * vars.chi)) * vars.h[i][j] * tr_A2 + (1.0 / (6 * vars.chi)) * vars.K * vars.A[i][j];
            FOR1(k)
            {
                Eij[i][j] -= (1.0 / (vars.chi)) * (vars.A[k][i] * d1.shift[m][j] - vars.A[k][j] * d1.shift[k][i]);
            }
        }
        return Eij;
    }

    //For computing the magnetic part of the Chern-Simons term B_{ij}.
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static Tensor<2, data_t>
    compute_magnetic_term(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        using namespace TensorAlgebra;

        Tensor<2, data_t> Bij = 0.0; //the Chern-Simons electric term

        const auto epsilon3_LUU = Weyl4::compute_epsilon3_LUU(vars, h_UU);
        const auto covd_Aij_LLL = compute_covd_Aij(vars, d1, h_UU, chris);

        //Finally compute the magnetic term of Chern-Simons term
        FOR2(i, j)
        {
            FOR2(s, m)
            {
                Bij[i][j] = epsilon3_LUU[i][s][m] * covd_Aij_LLL[j][m][s] + epsilon3_LUU[l][s][m] * covd_Aij_LLL[i][m][s];
            }
            
        }
        return Bij;
    }
    
    //For computing the Pontryagin density *RR, i.e. the Chern-Simons term.
    template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
    static modifiedscalar_t<data_t>
    compute_chern_simons(const vars_t<data_t> &vars,
        const vars_t<Tensor<1, data_t>> &d1, const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris)
    {
        using namespace TensorAlgebra;

        modifiedscalar_t<data_t> out;

        const auto E_ij = compute_chern_simons_electric_term(vars, d1, d2, h_UU, chris);
        const auto B_ij = compute_magnetic_term(vars, d1, d2, h_UU, chris);

        //Finally compute *RR
        FOR4(i, j, k, l)
        {
            out.starR_R = - 8.0 * vars.chi * vars.chi * h_UU[k][i] * h_UU[l][j] * B_ij[k][l] * E_ij[i][j];
        }
        return out;
    }

};

#endif /* CCZ4GEOMETRYMODIFIEDGR_HPP_ */
