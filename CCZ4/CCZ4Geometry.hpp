//This file calculates CCZ4 geometric quantities (or a similar 3+1 split.
//It assume input in terms of the struct CCZ4::vars_t defined in CCZ4.hpp
#ifndef CCZ4GEOMETRY_HPP_
#define CCZ4GEOMETRY_HPP_

#include "TensorAlgebra.hpp"

template <class data_t>
struct chris_t
{
        tensor<3, data_t> LLL; //3 lower indices
        tensor<3, data_t> ULL;
        tensor<1, data_t> contracted;  //chris contracted
};

template <class data_t>
struct ricci_t
{
     tensor<2, data_t> LL; //Ricci with two indices down
     data_t            scalar; //Ricci scalar
};

class CCZ4Geometry
{
    public:
    template <class data_t, template <typename> class vars_t>
    static tensor<2, data_t>
    compute_inverse_metric(const vars_t<data_t> &vars)
    {
        data_t deth = vars.h[0][0]*vars.h[1][1]*vars.h[2][2] + 2*vars.h[0][1]*vars.h[0][2]*vars.h[1][2] - vars.h[0][0]*vars.h[1][2]*vars.h[1][2] - vars.h[1][1]*vars.h[0][2]*vars.h[0][2] - vars.h[2][2]*vars.h[0][1]*vars.h[0][1];
        tensor<2, data_t> h_UU;
        {
                h_UU[0][0] = (vars.h[1][1]*vars.h[2][2] - vars.h[1][2]*vars.h[1][2]) / deth;
                h_UU[0][1] = (vars.h[0][2]*vars.h[1][2] - vars.h[0][1]*vars.h[2][2]) / deth;
                h_UU[0][2] = (vars.h[0][1]*vars.h[1][2] - vars.h[0][2]*vars.h[1][1]) / deth;
                h_UU[1][1] = (vars.h[0][0]*vars.h[2][2] - vars.h[0][2]*vars.h[0][2]) / deth;
                h_UU[1][2] = (vars.h[0][1]*vars.h[0][2] - vars.h[0][0]*vars.h[1][2]) / deth;
                h_UU[2][2] = (vars.h[0][0]*vars.h[1][1] - vars.h[0][1]*vars.h[0][1]) / deth;
                h_UU[1][0] = h_UU[0][1];
                h_UU[2][0] = h_UU[0][2];
                h_UU[2][1] = h_UU[1][2];
        }
        return h_UU;
    };

    template <class data_t, template <typename> class vars_t>
    static chris_t<data_t>
    compute_christoffel(
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const tensor<2, data_t>& h_UU
    )
    {
        chris_t<data_t> out;

        FOR3(i,j,k)
        {
            out.LLL[i][j][k] = 0.5*(d1[k].h[j][i] + d1[j].h[k][i] - d1[i].h[j][k]);
        }
        FOR3(i,j,k)
        {
            out.ULL[i][j][k] = 0;
            FOR1(l)
            {
                 out.ULL[i][j][k] += h_UU[i][l]*out.LLL[l][j][k];
            }
        }

        // Technically we can write out.contracted[i] += h_UU[j][k]*chris.ULL[i][j][k],
        // but sometimes people write:
        // out.contracted[i] += h_UU[i][j]*h_UU[k][l]*d1[l].h[k][j];
        // In theory h_UU[j][k]*d1[i].h[j][k] should be zero due to det h = 1
        // but in practice this term can deviate from zero.
        // For PRL 116, 071102 we used the former and it seemed to work well.
        FOR1(i)
        {
            out.contracted[i] = 0;
            FOR2(j,k)
            {
                 out.contracted[i] += h_UU[j][k]*out.ULL[i][j][k];
            }
            //FOR3(j,k,l)
            //{
            //     //out.contracted[i] += h_UU[i][j]*h_UU[k][l]*d1[l].h[k][j];
            //}
        }

        return out;
    };

    template <class data_t, template <typename> class vars_t>
    static ricci_t<data_t>
    compute_ricci_Z(
        const vars_t<data_t> &vars,
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
        const tensor<2, data_t>& h_UU,
        const chris_t<data_t>& chris,
        const tensor<1, data_t>& Z_over_chi
    )
    {
        ricci_t<data_t> out;

        data_t boxtildechi = 0;

        tensor<2, data_t> covdtilde2chi;
        FOR2(k,l)
        {
            covdtilde2chi[k][l] = d2[k][l].chi;
            FOR1(m)
            {
                 covdtilde2chi[k][l] -= chris.ULL[m][k][l]*d1[m].chi;
            }
        }

        FOR2(k,l)
        {
            boxtildechi += h_UU[k][l]*covdtilde2chi[k][l];
        }

        data_t dchi_dot_dchi = 0;
        {
            FOR2(m,n)
            {
                 dchi_dot_dchi += h_UU[m][n]*d1[m].chi*d1[n].chi;
            }
        }

        FOR2(i,j)
        {
            data_t ricci_tilde = 0;
            FOR1(k)
            {
                 // Trick: For CCZ4, we can add Z terms to ricci by changing Gamma to chrisvec
                 // This way of writing it allows the user to pass Z/chi = {0};
                 ricci_tilde += 0.5*(vars.h[k][i]*d1[j].Gamma[k] + vars.h[k][j]*d1[i].Gamma[k]);
                 ricci_tilde += 0.5*(vars.Gamma[k] - 2*Z_over_chi[k])*(chris.LLL[i][j][k] + chris.LLL[j][i][k]);
                 FOR1(l)
                 {
                        ricci_tilde -= 0.5*h_UU[k][l]*d2[k][l].h[i][j];
                        FOR1(m)
                        {
                             ricci_tilde += h_UU[l][m]*(chris.ULL[k][l][i]*chris.LLL[j][k][m] + chris.ULL[k][l][j]*chris.LLL[i][k][m] + chris.ULL[k][i][m]*chris.LLL[k][l][j]);
                        }
                 }
            }

            data_t ricci_chi = 0.5*((GR_SPACEDIM-2)*covdtilde2chi[i][j] + vars.h[i][j]*boxtildechi - ((GR_SPACEDIM-2)*d1[i].chi*d1[j].chi + GR_SPACEDIM*vars.h[i][j]*dchi_dot_dchi) / (2*vars.chi));

            data_t z_terms = 0;
            FOR1(k)
            {
                 z_terms += Z_over_chi[k]*(vars.h[i][k]*d1[j].chi + vars.h[j][k]*d1[i].chi - vars.h[i][j]*d1[k].chi + d1[k].h[i][j]*vars.chi);
            }

            out.LL[i][j] = (ricci_chi + vars.chi*ricci_tilde + z_terms) / vars.chi;
        }

        out.scalar = vars.chi*compute_trace(out.LL, h_UU);

        return out;
    }

    template <class data_t, template <typename> class vars_t>
    static ricci_t<data_t>
    compute_ricci(
        const vars_t<data_t> &vars,
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
        const tensor<2, data_t>& h_UU,
        const chris_t<data_t>& chris
    )
    {
         tensor<1,data_t> Z0;
         FOR1(i) Z0[i] = 0; //TODO: fix the array constructor in tensor
         return compute_ricci_Z(vars, d1, d2, h_UU, chris, Z0);
    }

};

#endif /* CCZ4GEOMETRY_HPP_ */
