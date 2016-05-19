/**
  * BOOSTED SCHWARZSCHILD BLACK HOLE
  * Baumgarte & Shapiro, pp. 73-74
  * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
  */

class BoostedBH
{

public:
	const Real mass;
 	const std::vector<Real> center;
 	const std::vector<Real> momentum;

	BoostedBH(Real m, std::vector<Real> const& center, std::vector<Real> const& momentum);

	// conformal factor
	Real psi(Real x, Real y, Real z) const;

	// extrinsic curvature
	void Aij(Real x, Real y, Real z, Real out[3][3]) const;

private:
	Real center_dist(Real x, Real y, Real z) const;
	Real psi0(Real r) const;
	Real psi2(Real r, Real cos_theta) const;
	Real psi2_0(Real r) const;
	Real psi2_2(Real r) const;

};

/* PUBLIC */

BoostedBH::BoostedBH(Real mass, std::vector<Real> const& center, std::vector<Real> const& momentum) :
	mass (mass),
	center (center),
	momentum (momentum)
{
	
}

Real
BoostedBH::psi(Real x, Real y, Real z) const
{
	const Real r = center_dist(x,y,z);
	const Real cos_theta = (z - center[2]) / r;
	const Real P_squared = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2];
	return psi0(r) + P_squared * psi2(r, cos_theta) / (mass * mass);
}

void
BoostedBH::Aij(Real x, Real y, Real z, Real out[3][3]) const
{
	const Real r = center_dist(x,y,z);
	const Real l[3] = { (x - center[0]) / r, (y - center[1]) / r, (z - center[2]) / r };
	const Real l_dot_p = l[0] * momentum[0] + l[1] * momentum[1] + l[2] * momentum[2];

	for (int i = 0; i < 3; ++i)
	{
		for (int j = i; j < 3; ++j)
		{
			const Real delta = (i == j) ? 1 : 0;
			out[i][j] = 1.5 * (momentum[i] * l[j] + momentum[j] * l[i] - (delta - l[i] * l[j]) * l_dot_p) / (r * r);
		}
	}
}

/* PRIVATE */

Real
BoostedBH::center_dist(Real x, Real y, Real z) const
{
	Real r = std::sqrt((x - center[0]) * (x - center[0]) + (y - center[1]) * (y - center[1]) + (z - center[2]) * (z - center[2]));
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
	return mass / (2 * r);
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

