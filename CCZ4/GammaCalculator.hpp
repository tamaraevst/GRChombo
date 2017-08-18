//Last update K Clough 23.02.2017

#ifndef GAMMACALCULATOR_HPP_
#define GAMMACALCULATOR_HPP_

#include "simd.hpp"
#include "VarsTools.hpp"
#include "GRInterval.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "CCZ4.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "Cell.hpp"
#include <array>

class GammaCalculator
{
//Use the variable definition in CCZ4
template<class data_t>
using Vars=CCZ4::Vars<data_t>;

protected:
    const FourthOrderDerivatives m_deriv;//!< An object for calculating derivatives of the variables
    double m_dx;

public:
    GammaCalculator(double a_dx) :  m_dx (a_dx), m_deriv (a_dx) {}

    template <class data_t>
    void compute(Cell<data_t> current_cell)
    {
        //copy data from chombo gridpoint into local variables, and calc 1st derivs
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        using namespace TensorAlgebra;
        const auto h_UU = compute_inverse_sym(vars.h);
        const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

        //assign values of Gamma^k = h_UU^ij * \tilde{Gamma}^k_ij in the output FArrayBox
        current_cell.store_vars(chris.contracted, GRInterval<c_Gamma1, c_Gamma3>());
    }

};

#endif /* GAMMACALCULATOR_HPP_ */
