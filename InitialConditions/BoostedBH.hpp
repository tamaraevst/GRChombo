#ifndef BOOSTEDBH_HPP_
#define BOOSTEDBH_HPP_
/**
  * BOOSTED SCHWARZSCHILD BLACK HOLE
  * Baumgarte & Shapiro, pp. 73-74
  * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
  */

#include <vector>
#include "tensor.hpp"

class BoostedBH
{

public:
    struct params_t
    {
        double mass;
        std::vector<double> center;
        std::vector<double> momentum;
    };

    const params_t m_params;

	BoostedBH(params_t a_params);

	// conformal factor
	double psi_minus_one(double x, double y, double z) const;

	// extrinsic curvature
	tensor<2, double> Aij(double x, double y, double z) const;

private:
	double center_dist(double x, double y, double z) const;
	double psi0(double r) const;
	double psi2(double r, double cos_theta) const;
	double psi2_0(double r) const;
	double psi2_2(double r) const;

};
#endif /*BOOSTEDBH_HPP_*/
