#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "ComplexPotential.hpp"
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <iostream>

int main()
{
    //initial vars (i.e. the ones at the centre)
    std::vector<double> central_vars{-0.071, 0.0, 0.1, 0.0};

    //initialise storage arrays
    std::vector<std::vector<double>> initial_var_arrays{};
    std::vector<double> rhos{};

    //initialise parameter struct
    Potential::params_t potential_params{1.0, 0.0};

    //initialise RHS Class and solution observer
    BosonStarRHS boson_star_rhs(potential_params);
    BosonStarSolutionObserver<std::vector, std::vector<double>>
        sol_observer(initial_var_arrays, rhos);

    //do integration
    try
    {
        boost::numeric::odeint::integrate(boson_star_rhs, central_vars, 0.0,
        20.0, 0.01, sol_observer);
    }
    catch (const char* exception)
    {
        std::cout << exception << " Integration stopped at rho = " <<
            rhos.back() << "\n";
    }

    std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "rho\t\t\tpsi\t\t\talpha\n";
    for(int i=0; i<rhos.size(); ++i)
    {
        std::cout << rhos[i] << "\t\t" << initial_var_arrays[i][2] << "\t\t" <<
        initial_var_arrays[i][0] << "\n";
    }
    return 0;
}
