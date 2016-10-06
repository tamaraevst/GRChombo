#ifndef RELAXATIONCHI_HPP_
#define RELAXATIONCHI_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4EMTensorSF.hpp"
#include "VarsBase.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

class RelaxationChi
{
public:
    template <class data_t>
    struct vars_t : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later

        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;
        data_t phi;
        data_t PiM;

        vars_t();
    };

protected:
    const double m_relaxspeed;
    const FABDriverBase& m_driver;
    const FourthOrderDerivatives m_deriv;

public:
    RelaxationChi(const FABDriverBase& driver, double dx, double relaxspeed);

    template <class data_t>
    void compute(int ix, int iy, int iz);

protected:
    template <class data_t>
    vars_t<data_t> rhs_equation(
        const vars_t<data_t> &vars,
        const vars_t< tensor<1,data_t> > &d1,
        const vars_t< tensor<2,data_t> > &d2,
        const vars_t<data_t> &advec
    );

};

#include "RelaxationChi.impl.hpp"

#endif /* RELAXATIONCHI_HPP_ */
