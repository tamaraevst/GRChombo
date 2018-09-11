#include "BosonStar.hpp" //for BosonStar::params_t struct
#include "ComplexPotential.hpp" //for Potential::params_t struct
#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarIntegrator.hpp" //for BosonStarIntegrator class
#include "BosonStarBinarySearch.hpp" //for BosonStarBinarySearch class
#include <vector> //for initial_data_t and initial_state_t
#include <iostream> //for writing to cout
#include <cmath> //for std::pow and M_PI

//some typedefs to make it easier to understand what we're using each type for
typedef std::vector<double> initial_state_t;
template<class T>
using initial_data_t = std::vector<T>;
typedef initial_data_t<double> initial_grid_t;

int main()
{
    //Set parameters
    const BosonStar::params_t params_BosonStar{0.6, 1.0e-14, 1.0e-14,
        0.125, 200.0, std::pow(2.0,-51)};
    const Potential::params_t params_potential{1.0, 4.0 * M_PI * 100.0};
    const double alpha_central_min{-2.06043};
    const double alpha_central_max{-2.06042};

    //Construct integrator object
    BosonStarIntegrator<initial_data_t, initial_state_t> integrator(
        params_BosonStar, params_potential);

    //first calculate the "min" solution
    integrator.doIntegration(alpha_central_min);
    auto sol_min = integrator.getSolution();

    //now calculate the "max" solution
    integrator.doIntegration(alpha_central_max);
    auto sol_max = integrator.getSolution();

    //Construct binary search object
    BosonStarBinarySearch<initial_data_t, initial_state_t> binary_search(
        params_BosonStar, params_potential, sol_min, sol_max);

    //do shooting
    binary_search.shoot();

    //get solution out and print frequency
    auto sol = binary_search.getShootedSolution();
    double frequency = sol.get_frequency();
    std::cout << "Star Frequency = " << frequency << "\n";

    //print solution
    /*
    std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "rho\t\t\tpsi\t\t\tPsi\t\t\talpha\t\t\tbeta\n";
    for(int i = 0; i < static_cast<int>(sol.get_grid().size()); ++i)
    {
        std::cout << sol.get_grid()[i] << "\t\t" << sol.get_psi()[i] << "\t\t"
        << sol.get_Psi()[i] << "\t\t" << sol.get_alpha()[i] << "\t\t"
        << sol.get_beta()[i] <<"\n";
    }
    */
    return 0;

}
