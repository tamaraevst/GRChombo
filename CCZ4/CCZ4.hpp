#ifndef CCZ4_HPP_
#define CCZ4_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "VarsBase.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

class CCZ4
{
public:

    enum
    {
     USE_CCZ4,
     USE_BSSN
    };

    template <class data_t>
    struct Vars : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later
        using VarsBase<data_t>::define_symmetric_enum_mapping;

        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;

        Vars();
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
    int m_formulation;
    const FABDriverBase& m_driver;
    const FourthOrderDerivatives m_deriv;

public:
    CCZ4(const FABDriverBase& driver, params_t params, double dx, double sigma, int formulation = USE_CCZ4, double cosmological_constant = 0);

    template <class data_t>
    void compute(int ix, int iy, int iz);

protected:
    template <class data_t, template<typename> class vars_t>
    vars_t<data_t> rhs_equation(
        const vars_t<data_t> &vars,
        const vars_t< tensor<1,data_t> > &d1,
        const vars_t< tensor<2,data_t> > &d2,
        const vars_t<data_t> &advec
    );
};

#include "CCZ4.impl.hpp"

#endif /* CCZ4_HPP_ */
