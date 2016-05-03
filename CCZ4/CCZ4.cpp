#include "CCZ4.hpp"

void CCZ4_exec(
	const BaseFab<double>& in,
	BaseFab<double>& out,
	Box box,
	CCZ4_params params,
	double dx,
	double sigma
)
{
	const int *in_lo = in.loVect();
	const int *in_hi = in.hiVect();
	const int stride_x = 1;
	const int stride_y = in_hi[0]-in_lo[0]+1;
	const int stride_z = (in_hi[1]-in_lo[1]+1)*stride_y;

	const int *out_lo = out.loVect();
	const int *out_hi = out.hiVect();
	const int out_stride_x = 1;
	const int out_stride_y = out_hi[0]-out_lo[0]+1;
	const int out_stride_z = (out_hi[1]-out_lo[1]+1)*out_stride_y;

	const int *box_lo = box.loVect();
	const int *box_hi = box.hiVect();

	// dataPtr in Chombo does CH_assert bound check
	// which we don't want to do in a loop
	const double *__restrict__ in_ptr[c_NUM];
	for (int i = 0; i < c_NUM; ++i)
	{
		in_ptr[i] = in.dataPtr(i);
	}

	double *__restrict__ out_ptr[c_NUM];
	for (int i = 0; i < c_NUM; ++i)
	{
		out_ptr[i] = out.dataPtr(i);
	}

	// Copying out a local struct of CCZ4 variables is the same as applying a trivial 1-wide stencil
	stencil<1> local;
	local.offset[0] = 0; local.weight[0] = 1.0;

	stencil<4> d1_stencil;
	d1_stencil.offset[0] = -2; d1_stencil.weight[0] = +8.33333333333333333333e-2 / dx;
	d1_stencil.offset[1] = -1; d1_stencil.weight[1] = -6.66666666666666666667e-1 / dx;
	d1_stencil.offset[2] = +1; d1_stencil.weight[2] = +6.66666666666666666667e-1 / dx;
	d1_stencil.offset[3] = +2; d1_stencil.weight[3] = -8.33333333333333333333e-2 / dx;

	stencil<5> d2_stencil;
	d2_stencil.offset[0] = -2; d2_stencil.weight[0] = -8.33333333333333333333e-2 / (dx*dx);
	d2_stencil.offset[1] = -1; d2_stencil.weight[1] = +1.33333333333333333333e+0 / (dx*dx);
	d2_stencil.offset[2] =  0; d2_stencil.weight[2] = -2.50000000000000000000e+0 / (dx*dx);
	d2_stencil.offset[3] = +1; d2_stencil.weight[3] = +1.33333333333333333333e+0 / (dx*dx);
	d2_stencil.offset[4] = +2; d2_stencil.weight[4] = -8.33333333333333333333e-2 / (dx*dx);

	stencil<5> d1_upwind;
	d1_upwind.offset[0] = -1; d1_upwind.weight[0] = -2.50000000000000000000e-1 / dx;
	d1_upwind.offset[1] =  0; d1_upwind.weight[1] = -8.33333333333333333333e-1 / dx;
	d1_upwind.offset[2] = +1; d1_upwind.weight[2] = +1.50000000000000000000e+0 / dx;
	d1_upwind.offset[3] = +2; d1_upwind.weight[3] = -5.00000000000000000000e-1 / dx;
	d1_upwind.offset[4] = +3; d1_upwind.weight[4] = +8.33333333333333333333e-2 / dx;

	stencil<7> dissipation;
	dissipation.offset[0] = -3; dissipation.weight[0] = +1.56250e-2 / dx;
	dissipation.offset[1] = -2; dissipation.weight[1] = -9.37500e-2 / dx;
	dissipation.offset[2] = -1; dissipation.weight[2] = +2.34375e-1 / dx;
	dissipation.offset[3] =  0; dissipation.weight[3] = -3.12500e-1 / dx;
	dissipation.offset[4] = +1; dissipation.weight[4] = +2.34375e-1 / dx;
	dissipation.offset[5] = +2; dissipation.weight[5] = -9.37500e-2 / dx;
	dissipation.offset[6] = +3; dissipation.weight[6] = +1.56250e-2 / dx;

#pragma omp parallel for default(shared) collapse(3)
	for (int z = box_lo[2]; z <= box_hi[2]; ++z)
	for (int y = box_lo[1]; y <= box_hi[1]; ++y)
	for (int x = box_lo[0]; x <= box_hi[0]; ++x)
	{
		const int idx = stride_z*(z-in_lo[2]) + stride_y*(y-in_lo[1]) + (x-in_lo[0]);
		
		CCZ4_vars<double> vars = {};
		CCZ4_apply(local, in_ptr, idx, 1, vars);
		
		CCZ4_vars<double> d1[CH_SPACEDIM] = {};
		CCZ4_apply(d1_stencil, in_ptr, idx, stride_x, d1[0]);
		CCZ4_apply(d1_stencil, in_ptr, idx, stride_y, d1[1]);
		CCZ4_apply(d1_stencil, in_ptr, idx, stride_z, d1[2]);

		CCZ4_vars<double> d2[CH_SPACEDIM][CH_SPACEDIM] = {};

		// Repeated derivatives
		CCZ4_apply(d2_stencil, in_ptr, idx, stride_x, d2[0][0]);
		CCZ4_apply(d2_stencil, in_ptr, idx, stride_y, d2[1][1]);
		CCZ4_apply(d2_stencil, in_ptr, idx, stride_z, d2[2][2]);

		// Mixed derivatives
		CCZ4_apply2(d1_stencil, d1_stencil, in_ptr, idx, stride_y, stride_x, d2[0][1]);
		CCZ4_apply2(d1_stencil, d1_stencil, in_ptr, idx, stride_z, stride_x, d2[0][2]);
		CCZ4_apply2(d1_stencil, d1_stencil, in_ptr, idx, stride_z, stride_y, d2[1][2]);

		d2[1][0] = d2[0][1];
		d2[2][0] = d2[0][2];
		d2[2][1] = d2[1][2];

		CCZ4_vars<double> advec = {};
		int sign_x = vars.shift[0] > 0 ? 1 : -1;
		int sign_y = vars.shift[1] > 0 ? 1 : -1;
		int sign_z = vars.shift[2] > 0 ? 1 : -1;
		CCZ4_apply(d1_upwind, in_ptr, idx, stride_x*sign_x, advec, sign_x*vars.shift[0]);
		CCZ4_apply(d1_upwind, in_ptr, idx, stride_y*sign_y, advec, sign_y*vars.shift[1]);
		CCZ4_apply(d1_upwind, in_ptr, idx, stride_z*sign_z, advec, sign_z*vars.shift[2]);

		CCZ4_vars<double> rhs;
		CCZ4_rhs(vars, d1, d2, advec, rhs, params);

		CCZ4_vars<double> dssp = {};
		CCZ4_apply(dissipation, in_ptr, idx, stride_x, dssp);
		CCZ4_apply(dissipation, in_ptr, idx, stride_y, dssp);
		CCZ4_apply(dissipation, in_ptr, idx, stride_z, dssp);

		// TODO: I really do not like this, but cannot think of a better way to do it yet...
		const int out_idx = out_stride_z*(z-out_lo[2]) + out_stride_y*(y-out_lo[1]) + (x-out_lo[0]);
		out_ptr[c_chi][out_idx]    = rhs.chi      + sigma * dssp.chi;
		out_ptr[c_h11][out_idx]    = rhs.h[0][0]  + sigma * dssp.h[0][0];
		out_ptr[c_h12][out_idx]    = rhs.h[0][1]  + sigma * dssp.h[0][1];
		out_ptr[c_h13][out_idx]    = rhs.h[0][2]  + sigma * dssp.h[0][2];
		out_ptr[c_h22][out_idx]    = rhs.h[1][1]  + sigma * dssp.h[1][1];
		out_ptr[c_h23][out_idx]    = rhs.h[1][2]  + sigma * dssp.h[1][2];
		out_ptr[c_h33][out_idx]    = rhs.h[2][2]  + sigma * dssp.h[2][2];
		out_ptr[c_K][out_idx]      = rhs.K        + sigma * dssp.K;
		out_ptr[c_A11][out_idx]    = rhs.A[0][0]  + sigma * dssp.A[0][0];
		out_ptr[c_A12][out_idx]    = rhs.A[0][1]  + sigma * dssp.A[0][1];
		out_ptr[c_A13][out_idx]    = rhs.A[0][2]  + sigma * dssp.A[0][2];
		out_ptr[c_A22][out_idx]    = rhs.A[1][1]  + sigma * dssp.A[1][1];
		out_ptr[c_A23][out_idx]    = rhs.A[1][2]  + sigma * dssp.A[1][2];
		out_ptr[c_A33][out_idx]    = rhs.A[2][2]  + sigma * dssp.A[2][2];
		out_ptr[c_Gamma1][out_idx] = rhs.Gamma[0] + sigma * dssp.Gamma[0];
		out_ptr[c_Gamma2][out_idx] = rhs.Gamma[1] + sigma * dssp.Gamma[1];
		out_ptr[c_Gamma3][out_idx] = rhs.Gamma[2] + sigma * dssp.Gamma[2];
		out_ptr[c_Theta][out_idx]  = rhs.Theta    + sigma * dssp.Theta;
		out_ptr[c_lapse][out_idx]  = rhs.lapse    + sigma * dssp.lapse;
		out_ptr[c_shift1][out_idx] = rhs.shift[0] + sigma * dssp.shift[0];
		out_ptr[c_shift2][out_idx] = rhs.shift[1] + sigma * dssp.shift[1];
		out_ptr[c_shift3][out_idx] = rhs.shift[2] + sigma * dssp.shift[2];
		out_ptr[c_B1][out_idx]     = rhs.B[0]     + sigma * dssp.B[0];
		out_ptr[c_B2][out_idx]     = rhs.B[1]     + sigma * dssp.B[1];
		out_ptr[c_B3][out_idx]     = rhs.B[2]     + sigma * dssp.B[2];
	}
}

template <class data_t, int stencil_size>
void CCZ4_apply(
	const stencil<stencil_size> &s,
	const data_t *__restrict__ *in_ptr,
	int idx,
	int stride,
	CCZ4_vars<data_t> &out,
	double multiplier
)
{
	for (int i = 0; i < stencil_size; ++i)
	{
		int idx_shifted = idx + s.offset[i]*stride;

		out.chi      += multiplier*s.weight[i]*in_ptr[c_chi][idx_shifted];
		out.h[0][0]  += multiplier*s.weight[i]*in_ptr[c_h11][idx_shifted];
		out.h[0][1]  += multiplier*s.weight[i]*in_ptr[c_h12][idx_shifted];
		out.h[0][2]  += multiplier*s.weight[i]*in_ptr[c_h13][idx_shifted];
		out.h[1][1]  += multiplier*s.weight[i]*in_ptr[c_h22][idx_shifted];
		out.h[1][2]  += multiplier*s.weight[i]*in_ptr[c_h23][idx_shifted];
		out.h[2][2]  += multiplier*s.weight[i]*in_ptr[c_h33][idx_shifted];

		out.K        += multiplier*s.weight[i]*in_ptr[c_K][idx_shifted];
		out.A[0][0]  += multiplier*s.weight[i]*in_ptr[c_A11][idx_shifted];
		out.A[0][1]  += multiplier*s.weight[i]*in_ptr[c_A12][idx_shifted];
		out.A[0][2]  += multiplier*s.weight[i]*in_ptr[c_A13][idx_shifted];
		out.A[1][1]  += multiplier*s.weight[i]*in_ptr[c_A22][idx_shifted];
		out.A[1][2]  += multiplier*s.weight[i]*in_ptr[c_A23][idx_shifted];
		out.A[2][2]  += multiplier*s.weight[i]*in_ptr[c_A33][idx_shifted];

		out.Gamma[0] += multiplier*s.weight[i]*in_ptr[c_Gamma1][idx_shifted];
		out.Gamma[1] += multiplier*s.weight[i]*in_ptr[c_Gamma2][idx_shifted];
		out.Gamma[2] += multiplier*s.weight[i]*in_ptr[c_Gamma3][idx_shifted];

		out.Theta    += multiplier*s.weight[i]*in_ptr[c_Theta][idx_shifted];

		out.lapse    += multiplier*s.weight[i]*in_ptr[c_lapse][idx_shifted];
		out.shift[0] += multiplier*s.weight[i]*in_ptr[c_shift1][idx_shifted];
		out.shift[1] += multiplier*s.weight[i]*in_ptr[c_shift2][idx_shifted];
		out.shift[2] += multiplier*s.weight[i]*in_ptr[c_shift3][idx_shifted];

		out.B[0]     += multiplier*s.weight[i]*in_ptr[c_B1][idx_shifted];
		out.B[1]     += multiplier*s.weight[i]*in_ptr[c_B2][idx_shifted];
		out.B[2]     += multiplier*s.weight[i]*in_ptr[c_B3][idx_shifted];
	}

	out.h[1][0] = out.h[0][1];
	out.h[2][0] = out.h[0][2];
	out.h[2][1] = out.h[1][2];

	out.A[1][0] = out.A[0][1];
	out.A[2][0] = out.A[0][2];
	out.A[2][1] = out.A[1][2];
}

template <class data_t, int stencil_size1, int stencil_size2>
void CCZ4_apply2(
	const stencil<stencil_size1> &s1,
	const stencil<stencil_size2> &s2,
	const data_t *__restrict__ *in_ptr,
	int idx,
	int stride1,
	int stride2,
	CCZ4_vars<data_t> &out
)
{
	for (int i = 0; i < stencil_size1; ++i)
	{
		CCZ4_apply(s2, in_ptr, idx+s1.offset[i]*stride1, stride2, out, s1.weight[i]);
	}
}

template <class data_t, bool covariantZ4>
void CCZ4_rhs(
	const CCZ4_vars<data_t> &vars,
	const CCZ4_vars<data_t> (&d1)[CH_SPACEDIM],
	const CCZ4_vars<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
	const CCZ4_vars<data_t> &advec,
	CCZ4_vars<data_t> &rhs,
	CCZ4_params &params
)
{
	const double chi_regularised = vars.chi;

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

	tensor<3, data_t> chris_LLL;
	{
		FOR3(i,j,k)
		{
			chris_LLL[i][j][k] = 0.5*(d1[k].h[j][i] + d1[j].h[k][i] - d1[i].h[j][k]);
		}
	}

	tensor<3, data_t> chris = {};
	{
		FOR3(i,j,k)
		{
			FOR1(l)
			{
				chris[i][j][k] += h_UU[i][l]*chris_LLL[l][j][k];
			}
		}
	}

	// Technically we can write	chrisvec[i] = h_UU[j][k]*chris[i][j][k],
	// but this is not numerically stable: h_UU[j][k]*d1[i].h[j][k] should be zero
	// but in practice can be > O(1).
	tensor<1, data_t> chrisvec = {};
	{
		FOR4(i,j,k,l)
		{
			chrisvec[i] += h_UU[i][j]*h_UU[k][l]*d1[l].h[k][j];
		}
	}

	
	tensor<1, data_t> Z_over_chi;
	tensor<1, data_t> Z;
	{
		FOR1(i)
		{
			Z_over_chi[i] = 0.5*(vars.Gamma[i] - chrisvec[i]);
			Z[i] = vars.chi*Z_over_chi[i];
		}
	}

	data_t divshift = 0;
	data_t Z_dot_d1lapse = 0;
	{
		FOR1(k)
		{
			divshift += d1[k].shift[k]; 
			Z_dot_d1lapse += Z[k]*d1[k].lapse;
		}
	}

	data_t dchi_dot_dchi = 0;
	data_t dlapse_dot_dchi = 0;
	{
		FOR2(m,n)
		{
			dchi_dot_dchi += h_UU[m][n]*d1[m].chi*d1[n].chi;
			dlapse_dot_dchi += h_UU[m][n]*d1[m].lapse*d1[n].chi;
		}
	}

	tensor<2, data_t> covdtilde2chi;
	tensor<2, data_t> covdtilde2lapse;
	tensor<2, data_t> covd2lapse;
	{
		FOR2(k,l)
		{
			covdtilde2chi[k][l] = d2[k][l].chi;
			covdtilde2lapse[k][l] = d2[k][l].lapse;
			FOR1(m)
			{
				covdtilde2chi[k][l] -= chris[m][k][l]*d1[m].chi;
				covdtilde2lapse[k][l] -= chris[m][k][l]*d1[m].lapse;
			}
			covd2lapse[k][l] = vars.chi*covdtilde2lapse[k][l] + 0.5*(d1[k].lapse*d1[l].chi + d1[k].chi*d1[l].lapse - vars.h[k][l]*dlapse_dot_dchi);
		}
	}

	double tr_covd2lapse = -(GR_SPACEDIM/2.0)*dlapse_dot_dchi;
	{
		FOR1(i)
		{
			tr_covd2lapse -= vars.chi*chrisvec[i]*d1[i].lapse;
			FOR1(j)
			{
				tr_covd2lapse += h_UU[i][j]*(vars.chi*d2[i][j].lapse + d1[i].lapse*d1[j].chi);
			}
		}
	}

	// Trick: For CCZ4, we can add Z terms to ricci by changing Gamma to chrisvec
	tensor<2, data_t> ricci;
	{
		data_t boxtildechi = 0;
		FOR2(k,l)
		{
			boxtildechi += h_UU[k][l]*covdtilde2chi[k][l];
		}

		FOR2(i,j)
		{		
			data_t ricci_tilde = 0;
			FOR1(k)
			{
				ricci_tilde += 0.5*(vars.h[k][i]*d1[j].Gamma[k] + vars.h[k][j]*d1[i].Gamma[k] + chrisvec[k]*(chris_LLL[i][j][k] + chris_LLL[j][i][k]));
				FOR1(l)
				{
					ricci_tilde -= 0.5*h_UU[k][l]*d2[k][l].h[i][j];
					FOR1(m)
					{
						ricci_tilde += h_UU[l][m]*(chris[k][l][i]*chris_LLL[j][k][m] + chris[k][l][j]*chris_LLL[i][k][m] + chris[k][i][m]*chris_LLL[k][l][j]);
					}
				}
			}

			data_t ricci_chi = 0.5*((GR_SPACEDIM-2)*covdtilde2chi[i][j] + vars.h[i][j]*boxtildechi - ((GR_SPACEDIM-2)*d1[i].chi*d1[j].chi + GR_SPACEDIM*vars.h[i][j]*dchi_dot_dchi) / (2*chi_regularised));

			data_t ricci_Z = 0;
			FOR1(k)
			{
				ricci_Z += Z_over_chi[k]*(vars.h[i][k]*d1[j].chi + vars.h[j][k]*d1[i].chi - vars.h[i][j]*d1[k].chi + d1[k].h[i][j]*vars.chi);
			}

			ricci[i][j] = (ricci_chi + vars.chi*ricci_tilde + ricci_Z) / vars.chi;
		}
	}

	tensor<2, data_t> A_UU = {};
	{
		FOR2(i,j)
		{
			FOR2(k,l)
			{
				A_UU[i][j] += h_UU[i][k]*h_UU[j][l]*vars.A[k][l];
			}
		}
	}

	data_t tr_ricci = 0;
	data_t tr_AA = 0;
	{
		FOR2(i,j)
		{
			tr_ricci += vars.chi*h_UU[i][j]*ricci[i][j];
			tr_AA += A_UU[i][j]*vars.A[i][j];
		}
	}

	{
		rhs.chi = advec.chi + (2.0/GR_SPACEDIM)*vars.chi*(vars.lapse*vars.K - divshift);
		FOR2(i,j)
		{
			rhs.h[i][j] = advec.h[i][j] - 2.0*vars.lapse*vars.A[i][j] - (2.0/GR_SPACEDIM)*vars.h[i][j]*divshift;
			FOR1(k)
			{
				rhs.h[i][j] += vars.h[k][i]*d1[j].shift[k] + vars.h[k][j]*d1[i].shift[k];
			}
		}
	}

	{
		tensor<2, data_t> Adot_TF_expr;
		FOR2(i,j)
		{
			Adot_TF_expr[i][j] = -covd2lapse[i][j] + vars.chi*vars.lapse*ricci[i][j];
		}

		data_t Adot_TF_trace = 0;
		FOR2(i,j)
		{
			Adot_TF_trace += h_UU[i][j]*Adot_TF_expr[i][j];
		}

		tensor<2, data_t> Adot_TF;
		FOR2(i,j)
		{
			Adot_TF[i][j] = Adot_TF_expr[i][j] - (1.0/GR_SPACEDIM)*Adot_TF_trace*vars.h[i][j];
		}

		FOR2(i,j)
		{
			rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] + vars.A[i][j]*(vars.lapse*(vars.K - 2*vars.Theta) - (2.0/GR_SPACEDIM)*divshift);
			FOR1(k)
			{
				rhs.A[i][j] += vars.A[k][i]*d1[j].shift[k] + vars.A[k][j]*d1[i].shift[k];
				FOR1(l)
				{
					rhs.A[i][j] -= 2*vars.lapse*h_UU[k][l]*vars.A[i][k]*vars.A[l][j];
				}
			}
		}
	}

	data_t kappa1_lapse = params.kappa1 * (covariantZ4 ? 1.0 : vars.lapse);

	{
		data_t Thetadot = 0.5*vars.lapse*(tr_ricci - tr_AA + ((GR_SPACEDIM-1)/(double) GR_SPACEDIM)*vars.K*vars.K - 2*vars.Theta*vars.K) - 0.5*vars.Theta*kappa1_lapse*((GR_SPACEDIM+1) + params.kappa2*(GR_SPACEDIM-1)) - Z_dot_d1lapse;
		rhs.Theta = advec.Theta + Thetadot;

		rhs.K = advec.K + 2*Thetadot + vars.lapse*(tr_AA + (1.0/GR_SPACEDIM)*vars.K*vars.K) + kappa1_lapse*(1-params.kappa2)*vars.Theta + 2*Z_dot_d1lapse - tr_covd2lapse;
	}

	tensor<1, data_t> Gammadot;
	{
		FOR1(i)
		{
			Gammadot[i] = (2.0/GR_SPACEDIM)*(divshift*(chrisvec[i] + 2*params.kappa3*Z_over_chi[i]) - 2*vars.lapse*vars.K*Z_over_chi[i]) - 2*kappa1_lapse*Z_over_chi[i];
			FOR1(j)
			{
				Gammadot[i] += 2*h_UU[i][j]*(vars.lapse*d1[j].Theta - vars.Theta*d1[j].lapse)
					- 2*A_UU[i][j]*d1[j].lapse
					- vars.lapse*((2*(GR_SPACEDIM-1)/(double) GR_SPACEDIM)*h_UU[i][j]*d1[j].K + GR_SPACEDIM*A_UU[i][j]*d1[j].chi/chi_regularised)
					- (chrisvec[j] + 2*params.kappa3*Z_over_chi[j])*d1[j].shift[i];

				FOR1(k)
				{
					Gammadot[i] += 2*vars.lapse*chris[i][j][k]*A_UU[j][k]
						+ h_UU[j][k]*d2[j][k].shift[i]
						+ ((GR_SPACEDIM-2)/(double) GR_SPACEDIM)*h_UU[i][j]*d2[j][k].shift[k];
				}
			}
		}

		FOR1(i)
		{
			rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
		}
	}
	
	double eta = 1;

	{
		rhs.lapse = params.lapse_advec_coeff*advec.lapse - 2*vars.lapse*(vars.K - 2*vars.Theta);
		FOR1(i)
		{
			rhs.shift[i] = params.shift_advec_coeff*advec.shift[i] + params.shift_gamma_coeff*vars.B[i];
			rhs.B[i] = params.shift_advec_coeff*advec.B[i] + (1 - params.shift_advec_coeff)*advec.Gamma[i] + Gammadot[i] - params.beta_driver*eta*vars.B[i];
		}
	}
}
