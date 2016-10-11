#include "BoostedBH.hpp"
#include <cmath>

BoostedBH::BoostedBH(params_t a_params) : m_params (a_params){}

double
BoostedBH::psi_minus_one(double x, double y, double z) const
{
	const double r = center_dist(x,y,z);
	const double cos_theta = (z - m_params.center[2]) / r;
	const double P_squared = pow(m_params.momentum[0],2) + pow(m_params.momentum[1],2) + pow(m_params.momentum[2],2);
	return psi0(r) + P_squared * psi2(r, cos_theta) / (m_params.mass * m_params.mass);
}

tensor<2,double>
BoostedBH::Aij(double x, double y, double z) const
{
	const double r = center_dist(x,y,z);
	const double l[3] = { (x - m_params.center[0]) / r, (y - m_params.center[1]) / r, (z - m_params.center[2]) / r };
	const double l_dot_p = l[0] * m_params.momentum[0] + l[1] * m_params.momentum[1] + l[2] * m_params.momentum[2];

    tensor<2, double> out;

    FOR2(i,j)
    {
			const double delta = (i == j) ? 1 : 0;
			out[i][j] = 1.5 * (m_params.momentum[i] * l[j] + m_params.momentum[j] * l[i] - (delta - l[i] * l[j]) * l_dot_p) / (r * r);
	}
    return out;
}

/* PRIVATE */

double
BoostedBH::center_dist(double x, double y, double z) const
{
	double r = std::sqrt(pow(x - m_params.center[0],2) + pow(y - m_params.center[1],2) + pow(z - m_params.center[2],2));
	if (std::fabs(r) < 1e-6)
	{
		return 1e-6;
	}
	else
	{
		return r;
	}
}

double
BoostedBH::psi0(double r) const
{
	return m_params.mass / (2 * r);
}

double
BoostedBH::psi2(double r, double cos_theta) const
{
	return psi2_0(r) + psi2_2(r) * (1.5 * cos_theta * cos_theta - 0.5);
}

double
BoostedBH::psi2_0(double r) const
{
	const double F = psi0(r);
	const double FF = F * F;
	return std::pow(1 + F, -5) * (F / 8) * (FF * FF + 5 * F * FF + 10 * FF + 10 * F + 5);
}

double
BoostedBH::psi2_2(double r) const
{
	const double F = psi0(r);
	const double FF = F * F;
	return 0.05 * std::pow(1 + F, -5) * FF * (84 * F * FF * FF + 378 * FF * FF + 658 * F * FF
		+ 539 * FF + 192 * F + 15) + 4.2 * F * FF * std::log(F / (1 + F));
}
