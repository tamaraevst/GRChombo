#include "BosonStar.hpp"
#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "BosonStarSolution.hpp"
#include "BosonStarBinarySearch.hpp"
#include "ComplexPotential.hpp"
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <iostream>
#include <stdexcept>

//some typedefs to make it easier to understand what we're using each type for
typedef std::vector<double> initial_state_t;
template<class T>
using initial_data_t = std::vector<T>;
typedef initial_data_t<double> initial_grid_t;

int main()
{
    //Set parameters
    const BosonStar::params_t params_CSF{0.1};
    const Potential::params_t params_potential{1.0, 0.0};

    //identify all the BCs
    double alpha_central_min{-0.071};
    double alpha_central_max{-0.07};
    const double beta_central{0.0};
    const double Psi_central{0.0};

    //Set central BCs
    initial_state_t central_vars_min{alpha_central_min, beta_central,
        params_CSF.central_amplitude_CSF, Psi_central};
    initial_state_t central_vars_max{alpha_central_max, beta_central,
        params_CSF.central_amplitude_CSF, Psi_central};

    //Set integration error tolerances
    const double abs_error{1.0e-14};
    const double rel_error{1.0e-14};

    //Set initial step size and max radius
    const double drho{2.0e-7};
    const double max_radius{200.0};

    //initialise storage arrays
    initial_data_t<initial_state_t> initial_var_arrays_min{};
    initial_data_t<initial_state_t> initial_var_arrays_max{};
    initial_grid_t initial_grid_min{};
    initial_grid_t initial_grid_max{};

    //initialise RHS Class and solution observer. Since odeint makes a copy of
    //sol_observer, it should be safe to just instantiate it once here for all
    //the iterations in the while loop.
    BosonStarRHS boson_star_rhs(params_potential);
    BosonStarSolutionObserver<initial_data_t, initial_state_t>
        sol_observer_min(initial_var_arrays_min, initial_grid_min);
    BosonStarSolutionObserver<initial_data_t, initial_state_t>
        sol_observer_max(initial_var_arrays_max, initial_grid_max);

    using namespace boost::numeric::odeint;
    typedef runge_kutta_cash_karp54<initial_state_t> error_stepper_t;
    //do integration for min case
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>
            (abs_error, rel_error), boson_star_rhs, central_vars_min, 0.0,
            max_radius, drho, sol_observer_min);
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << " max radius = " <<
            initial_grid_min.back() << "\n";
    }

    //do integration for max case
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>
            (abs_error, rel_error), boson_star_rhs, central_vars_max, 0.0,
            max_radius, drho, sol_observer_max);
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << " max radius = " <<
            initial_grid_max.back() << "\n";
    }

    //Put arrays into solution objects
    BosonStarSolution<initial_data_t, initial_state_t>
        sol_min(initial_var_arrays_min, initial_grid_min);
    BosonStarSolution<initial_data_t, initial_state_t>
        sol_max(initial_var_arrays_max, initial_grid_max);

    //Construct BinarySearch object for shooting
    BosonStarBinarySearch<initial_data_t, initial_state_t> binary_search(
        params_CSF, params_potential, sol_min, sol_max);

    //do shooting
    binary_search.shoot(1e-15, max_radius);

    //get the solution out
    auto sol = binary_search.getSolution();

    std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "rho\t\t\tpsi\t\t\tPsi\t\t\talpha\t\t\tbeta\n";
    for(int i = 0; i < sol.get_grid().size(); ++i)
    {
        std::cout << sol.get_grid()[i] << "\t\t" << sol.get_psi()[i] << "\t\t"
        << sol.get_Psi()[i] << "\t\t" << sol.get_alpha()[i] << "\t\t"
        << sol.get_beta()[i] <<"\n";
    }
    return 0;
}
