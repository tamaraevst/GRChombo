#ifndef BOOSTEDBH_HPP_
#define BOOSTEDBH_HPP_
/**
  * BOOSTED SCHWARZSCHILD BLACK HOLE
  * Baumgarte & Shapiro, pp. 73-74
  * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
  */

class BoostedBH
{

public:
    struct params_t
    {
        Real mass;
        std::vector<Real> center;
        std::vector<Real> momentum;
    };

    const params_t m_params;

	BoostedBH(params_t a_params);

	// conformal factor
	Real psi_minus_one(Real x, Real y, Real z) const;

	// extrinsic curvature
	tensor<2, double> Aij(Real x, Real y, Real z) const;

private:
	Real center_dist(Real x, Real y, Real z) const;
	Real psi0(Real r) const;
	Real psi2(Real r, Real cos_theta) const;
	Real psi2_0(Real r) const;
	Real psi2_2(Real r) const;

};

/* PUBLIC */

BoostedBH::BoostedBH(params_t a_params) : m_params (a_params){}

Real
BoostedBH::psi_minus_one(Real x, Real y, Real z) const
{
	const Real r = center_dist(x,y,z);
	const Real cos_theta = (z - m_params.center[2]) / r;
	const Real P_squared = pow(m_params.momentum[0],2) + pow(m_params.momentum[1],2) + pow(m_params.momentum[2],2);
	return psi0(r) + P_squared * psi2(r, cos_theta) / (m_params.mass * m_params.mass);
}

tensor<2,double>
BoostedBH::Aij(Real x, Real y, Real z) const
{
	const Real r = center_dist(x,y,z);
	const Real l[3] = { (x - m_params.center[0]) / r, (y - m_params.center[1]) / r, (z - m_params.center[2]) / r };
	const Real l_dot_p = l[0] * m_params.momentum[0] + l[1] * m_params.momentum[1] + l[2] * m_params.momentum[2];

    tensor<2, double> out;

    FOR2(i,j)
    {
			const Real delta = (i == j) ? 1 : 0;
			out[i][j] = 1.5 * (m_params.momentum[i] * l[j] + m_params.momentum[j] * l[i] - (delta - l[i] * l[j]) * l_dot_p) / (r * r);
	}
    return out;
}

/* PRIVATE */

Real
BoostedBH::center_dist(Real x, Real y, Real z) const
{
	Real r = std::sqrt(pow(x - m_params.center[0],2) + pow(y - m_params.center[1],2) + pow(z - m_params.center[2],2));
	if (std::fabs(r) < 1e-6)
	{
		return 1e-6;
	}
	else
	{
		return r;
	}
}

Real
BoostedBH::psi0(Real r) const
{
	return m_params.mass / (2 * r);
}

Real
BoostedBH::psi2(Real r, Real cos_theta) const
{
	return psi2_0(r) + psi2_2(r) * (1.5 * cos_theta * cos_theta - 0.5);
}

Real
BoostedBH::psi2_0(Real r) const
{
	const Real F = psi0(r);
	const Real FF = F * F;
	return std::pow(1 + F, -5) * (F / 8) * (FF * FF + 5 * F * FF + 10 * FF + 10 * F + 5);
}

Real
BoostedBH::psi2_2(Real r) const
{
	const Real F = psi0(r);
	const Real FF = F * F;
	return 0.05 * std::pow(1 + F, -5) * FF * (84 * F * FF * FF + 378 * FF * FF + 658 * F * FF
		+ 539 * FF + 192 * F + 15) + 4.2 * F * FF * std::log(F / (1 + F));
}
#endif /*BOOSTEDBH_HPP_*/
