#include "BosonStar.hpp" //for BosonStar::params_t struct
#include "ComplexPotential.hpp" //for Potential::params_t struct
#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarIntegrator.hpp" //for BosonStarIntegrator class
#include "BosonStarBinarySearch.hpp" //for BosonStarBinarySearch class
#include "BosonStarIsotropicSolution.hpp" //for BosonStarIsotropicSolution class
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
    const BosonStar::params_t params_BosonStar{0.3, 1.0e-14, 1.0e-14,
        std::pow(2.0,-6), 200.0, std::pow(2.0,-51)};
    const Potential::params_t params_potential{1.0, 4.0 * M_PI * 50.0};

    try
    {
        BosonStarBinarySearch<initial_data_t, initial_state_t> binary_search(
        params_BosonStar, params_potential);

        binary_search.shoot();
        auto sol = binary_search.getShootedSolution();

        BosonStarIsotropicSolution<initial_data_t, initial_state_t>
            isotropic_sol(params_BosonStar, params_potential);

        isotropic_sol.makeFromPolarArealSolution(sol, 200.0);

        std::cout.precision(16);
        std::cout << std::fixed;
        std::cout << "R\t\t\tphi\t\t\tchi\t\t\tlapse" << "\n";
        double R;
        for(int i = 0; i < 400; ++i)
        {
            R = i/2.0;
            std::cout << R << "\t" << isotropic_sol.m_phi(R) << "\t" <<
            isotropic_sol.m_chi(R) << "\t" << isotropic_sol.m_lapse(R) << "\n";
        }

        /*
        //print solution
        std::cout.precision(10);
        std::cout << std::fixed;
        std::cout << "rho\t\t\tpsi\t\t\tPsi\t\t\talpha\t\t\tbeta\n";
        for(int i = 0; i < static_cast<int>(sol.get_grid().size()); ++i)
        {
            std::cout << sol.get_grid()[i] << "\t\t" << sol.get_psi()[i] << "\t\t"
            << sol.get_Psi()[i] << "\t\t" << sol.get_alpha()[i] << "\t\t"
            << sol.get_beta()[i] << "\n";
        }
        */
    }
    catch (std::exception &exception)
    {
        std::cout << exception.what() << "\n";
    }

    return 0;

}
