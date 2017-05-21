//Last update K Clough 23.02.2017

#ifndef GAMMACALCULATOR_HPP_
#define GAMMACALCULATOR_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
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
    void compute(Cell current_cell)
    {
        //copy data from chombo gridpoint into local variables
        Vars<data_t> vars;
        current_cell.local_vars(vars);

        //work out first derivatives of variables on grid
        Vars< tensor<1, data_t> > d1;
        FOR1(idir) m_deriv.diff1(d1, current_cell, idir);

        using namespace TensorAlgebra;
        auto h_UU = compute_inverse(vars.h);
        auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

        //assign values of Gamma^k = h_UU^ij * \tilde{Gamma}^k_ij
        FOR1(i)
        {
            vars.Gamma[i] = 0.0;
            FOR2(j,k)
            {
                vars.Gamma[i] += h_UU[j][k] * chris.ULL[i][j][k];
            }
        }

        //Write the rhs into the output FArrayBox
        current_cell.store_vars(vars.Gamma[0], c_Gamma1);
        current_cell.store_vars(vars.Gamma[1], c_Gamma2);
        current_cell.store_vars(vars.Gamma[2], c_Gamma3);
    }

};

#endif /* GAMMACALCULATOR_HPP_ */
