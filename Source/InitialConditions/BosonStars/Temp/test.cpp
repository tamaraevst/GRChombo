#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
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
    //initial vars (i.e. the ones at the centre)
    initial_state_t central_vars{-0.159, 0.0, 0.1, 0.0};

    //initialise storage arrays
    initial_data_t<initial_state_t> initial_var_arrays{};
    initial_grid_t rhos{};

    //initialise parameter struct
    Potential::params_t potential_params{1.0, 0.0};

    //initialise RHS Class and solution observer
    BosonStarRHS boson_star_rhs(potential_params);
    int num_psi_roots{0};
    BosonStarSolutionObserver<initial_data_t, initial_state_t>
        sol_observer(initial_var_arrays, rhos, num_psi_roots);

    using namespace boost::numeric::odeint;
    typedef runge_kutta_cash_karp54<initial_state_t> error_stepper_t;
    //do integration
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>(1.0e-14, 1.0e-12),
            boson_star_rhs, central_vars, 0.0, 100.0, 0.01, sol_observer);
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << " rho_max = " <<
            rhos.back() << "\n";
    }
    std::cout << "The number of roots in psi is " << num_psi_roots << "\n";
    /*std::cout.precision(10);
    std::cout << std::fixed;
    std::cout << "rho\t\t\tpsi\t\t\tPsi\t\t\talpha\t\t\tbeta\n";
    for(int i=0; i<rhos.size(); ++i)
    {
        std::cout << rhos[i] << "\t\t" << initial_var_arrays[i][2] << "\t\t" <<
        initial_var_arrays[i][3] << "\t\t" << initial_var_arrays[i][0] << "\t\t"
        << initial_var_arrays[i][1] <<"\n";
    }*/
    return 0;
}
