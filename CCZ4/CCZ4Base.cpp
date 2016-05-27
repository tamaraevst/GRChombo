#ifndef CCZ4BASE_HPP_
#define CCZ4BASE_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"

class CCZ4Base
{
protected:
    template <class data_t>
    struct vars_t
    {
        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;
    };

public:
    template <class data_t>
    struct inverse_metric_t
    {
        data_t deth;
        tensor<2, data_t> h_UU;
    };

    template <class data_t>
    struct christoffel_t
    {
        tensor<3, data_t> chris_LLL;
        tensor<3, data_t> chris;
        tensor<1, data_t> chrisvec;
    };

    template <class data_t>
    ALWAYS_INLINE
    inverse_metric_t<data_t> compute_inverse_metric(const vars_t<data_t> &vars)
    {
        basic_geometry_t<data_t> out;

        out.deth = vars.h[0][0]*vars.h[1][1]*vars.h[2][2] + 2*vars.h[0][1]*vars.h[0][2]*vars.h[1][2] - vars.h[0][0]*vars.h[1][2]*vars.h[1][2] - vars.h[1][1]*vars.h[0][2]*vars.h[0][2] - vars.h[2][2]*vars.h[0][1]*vars.h[0][1];
        
        out.h_UU[0][0] = (vars.h[1][1]*vars.h[2][2] - vars.h[1][2]*vars.h[1][2]) / out.deth;
        out.h_UU[0][1] = (vars.h[0][2]*vars.h[1][2] - vars.h[0][1]*vars.h[2][2]) / out.deth;
        out.h_UU[0][2] = (vars.h[0][1]*vars.h[1][2] - vars.h[0][2]*vars.h[1][1]) / out.deth;
        out.h_UU[1][1] = (vars.h[0][0]*vars.h[2][2] - vars.h[0][2]*vars.h[0][2]) / out.deth;
        out.h_UU[1][2] = (vars.h[0][1]*vars.h[0][2] - vars.h[0][0]*vars.h[1][2]) / out.deth;
        out.h_UU[2][2] = (vars.h[0][0]*vars.h[1][1] - vars.h[0][1]*vars.h[0][1]) / out.deth;
        out.h_UU[1][0] = h_UU[0][1];
        out.h_UU[2][0] = h_UU[0][2];
        out.h_UU[2][1] = h_UU[1][2];
        

        return out;
    };

    template <class data_t>
    ALWAYS_INLINE
    christoffel_t<data_t> compute_christoffel(
        const vars_t<data_t> &vars,
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const inverse_metric_t<data_t> inv
    )
    {
        christoffel_t<data_t> out;

        FOR3(i,j,k)
        {
            out.chris_LLL[i][j][k] = 0.5*(d1[k].h[j][i] + d1[j].h[k][i] - d1[i].h[j][k]);
        }

        FOR3(i,j,k)
        {
            out.chris[i][j][k] = 0;
            FOR1(l)
            {
                out.chris[i][j][k] += inv.h_UU[i][l]*out.chris_LLL[l][j][k];
            }
        }

        // Technically we can write chrisvec[i] = h_UU[j][k]*chris[i][j][k],
        // but this is not numerically stable: h_UU[j][k]*d1[i].h[j][k] should be zero
        // but in practice can be > O(1).
        FOR1(i)
        {
            out.chrisvec[i] = 0;
            FOR3(j,k,l)
            {
                out.chrisvec[i] += inv.h_UU[i][j]*inv.h_UU[k][l]*d1[l].h[k][j];
            }
        }

        return out;
    }

    struct params_t
    {
        double kappa1;
        double kappa2;
        double kappa3;
        double shift_gamma_coeff;
        double lapse_advec_coeff;
        double shift_advec_coeff;
        double beta_driver;
    };

protected:
    params_t m_params;
    double m_dx;
    double m_sigma;

    const double *m_in_ptr[c_NUM];
    const int *m_in_lo;
    const int *m_in_hi;
    int m_stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];
    
public:
    CCZ4(params_t params, double dx, double sigma);
    void execute(const FArrayBox& in, FArrayBox& out);

protected:
    template <class data_t>
    void compute(int x, int y, int z);

    template <class data_t>
    void demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out);

    template <class data_t>
    void local_vars(int idx, vars_t<data_t>& out);

    template <class data_t>
    void diff1(int idx, int stride, vars_t<data_t>& out);

    template <class data_t>
    void diff2(int idx, int stride, vars_t<data_t>& out);

    template <class data_t>
    void mixed_diff2(int idx, int stride1, int stride2, vars_t<data_t>& out);

    template <class data_t>
    void advection(int idx, const tensor<1, data_t>& shift, vars_t<data_t>& out);

    template <class data_t>
    void dissipation(int idx, vars_t<data_t>& out);

    template <class data_t>
    void rhs_equation(const vars_t<data_t> &vars,
             const vars_t<data_t> (&d1)[CH_SPACEDIM],
             const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
             const vars_t<data_t> &advec,
             vars_t<data_t> &rhs);
};