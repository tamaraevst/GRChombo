#ifndef _CCZ4_HPP_
#define _CCZ4_HPP_

#include "GRUtils.H"
#include "BaseFab.H"

enum
{
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

template <class data_t>
struct CCZ4_vars
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

struct CCZ4_params
{
	double kappa1;
	double kappa2;
	double kappa3;
	double shift_gamma_coeff;
	double lapse_advec_coeff;
	double shift_advec_coeff;
	double beta_driver;
};

// Main RHS driver, loops through a box, apply stencils, etc.
// Can make this explicitly SIMD if required
void CCZ4_exec(
	const BaseFab<double>& in,
	BaseFab<double>& out,
	Box box,
	CCZ4_params params,
	double dx,
	double sigma
);

// General stencil applicator
template <class data_t, int stencil_size>
void CCZ4_apply(
	const stencil<stencil_size> &s,
	const data_t *__restrict__ *in_ptr,
	int idx,
	int stride,
	CCZ4_vars<data_t> &out,
	double multiplier = 1
);

// Apply tensor product of stencils, e.g. for mixed second derivatives
template <class data_t, int stencil_size1, int stencil_size2>
void CCZ4_apply2(
	const stencil<stencil_size1> &s1,
	const stencil<stencil_size2> &s2,
	const data_t *__restrict__ *in_ptr,
	int idx,
	int stride1,
	int stride2,
	CCZ4_vars<data_t> &out
);

// Actual RHS equations
template <class data_t, bool covariantZ4 = true>
void CCZ4_rhs(
	const CCZ4_vars<data_t> &vars,
	const CCZ4_vars<data_t> (&d1)[CH_SPACEDIM],
	const CCZ4_vars<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
	const CCZ4_vars<data_t> &advec,
	CCZ4_vars<data_t> &rhs,
	CCZ4_params &params
);

#endif /* _CCZ4_HPP_ */
