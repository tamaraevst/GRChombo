#if !defined(BOOSTEDBH_HPP_)
#error "This file should only be included through CCZ4.hpp"
#endif

#ifndef BOOSTEDBH_IMPL_HPP_
#define BOOSTEDBH_IMPL_HPP_

#include "BoostedBH.hpp"
#include <cmath>

inline
BoostedBH::BoostedBH(params_t a_params) : m_params (a_params){}

template <class data_t>
data_t
BoostedBH::psi_minus_one(Coordinates<data_t> coords) const
{
	const data_t r = center_dist(coords);
	const data_t cos_theta = (coords.z - m_params.center[2]) / r;
	const data_t P_squared = pow(m_params.momentum[0],2) + pow(m_params.momentum[1],2) + pow(m_params.momentum[2],2);
	return psi0(r) + P_squared * psi2(r, cos_theta) / (m_params.mass * m_params.mass);
}

template <class data_t>
tensor<2,data_t>
BoostedBH::Aij(Coordinates<data_t> coords) const
{
	const data_t r = center_dist(coords);
	const data_t l[3] = { (coords.x - m_params.center[0]) / r, (coords.y - m_params.center[1]) / r, (coords.z - m_params.center[2]) / r };
	const data_t l_dot_p = l[0] * m_params.momentum[0] + l[1] * m_params.momentum[1] + l[2] * m_params.momentum[2];

    tensor<2, data_t> out;

    FOR2(i,j)
    {
			const double delta = (i == j) ? 1 : 0;
			out[i][j] = 1.5 * (m_params.momentum[i] * l[j] + m_params.momentum[j] * l[i] - (delta - l[i] * l[j]) * l_dot_p) / (r * r);
	}
    return out;
}

/* PRIVATE */

template <class data_t>
data_t
BoostedBH::center_dist(Coordinates<data_t> coords) const
{
	data_t r = sqrt(  pow(coords.x - m_params.center[0],(decltype(coords.x))2)
                    + pow(coords.y - m_params.center[1],2)
                    + pow(coords.z - m_params.center[2],2));

    double minimum_r = 1e-6;
    auto r_is_too_small = simd_compare_lt(r, minimum_r);
    return simd_conditional(r_is_too_small, minimum_r, r);
}

template <class data_t>
data_t
BoostedBH::psi0(data_t r) const
{
	return m_params.mass / (2 * r);
}

template <class data_t>
data_t
BoostedBH::psi2(data_t r, data_t cos_theta) const
{
	return psi2_0(r) + psi2_2(r) * (1.5 * cos_theta * cos_theta - 0.5);
}

template <class data_t>
data_t
BoostedBH::psi2_0(data_t r) const
{
	const data_t F = psi0(r);
	const data_t FF = F * F;
	return pow(1 + F, (decltype(F))-5) * (F / 8) * (FF * FF + 5 * F * FF + 10 * FF + 10 * F + 5);
}

template <class data_t>
data_t
BoostedBH::psi2_2(data_t r) const
{
	const data_t F = psi0(r);
	const data_t FF = F * F;
	return 0.05 * pow(1 + F, (decltype(F))-5) * FF * (84 * F * FF * FF + 378 * FF * FF + 658 * F * FF
		+ 539 * FF + 192 * F + 15) + 4.2 * F * FF * log(F / (1 + F));
}

#endif /* BOOSTEDBH_IMPL_HPP_ */
