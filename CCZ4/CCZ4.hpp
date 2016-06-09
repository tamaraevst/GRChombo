#ifndef CCZ4_HPP_
#define CCZ4_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"

#include "user_enum.hpp" //This files needs c_NUM - total number of components

#include <array>

class CCZ4
{
public:
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

        vars_t(){};

        template <class arr_t>
        vars_t(const arr_t& in);
    };

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
    const params_t m_params;
    const double m_sigma;
    double m_cosmological_constant;
    const FABDriverBase& m_driver;
    const FourthOrderDerivatives m_deriv;

public:
    CCZ4(const FABDriverBase& driver, params_t params, double dx, double sigma, double cosmological_constant = 0);

    template <class data_t>
    void compute(int ix, int iy, int iz);

protected:
    template <class data_t>
    vars_t<data_t> rhs_equation(
        const vars_t<data_t> &vars,
        const vars_t<data_t> (&d1)[CH_SPACEDIM],
        const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
        const vars_t<data_t> &advec
    );
};

#include "CCZ4.impl.hpp"

#endif /* CCZ4_HPP_ */
