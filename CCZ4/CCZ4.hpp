#ifndef CCZ4_HPP_
#define CCZ4_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "GRUtils.H"

enum {
    c_chi,

    c_h,
    c_h11 = c_h,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A,
    c_A11 = c_A,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma,
    c_Gamma1 = c_Gamma,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift,
    c_shift1 = c_shift,
    c_shift2,
    c_shift3,

    c_B,
    c_B1 = c_B,
    c_B2,
    c_B3,

    c_NUM
};

class CCZ4
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

#endif /* CCZ4_HPP_ */
