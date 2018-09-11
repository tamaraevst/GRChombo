/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

#include <cmath> //for std::exp and std::abs
#include <iostream> //TODO: remove after debugging
#include <boost/math/interpolators/barycentric_rational.hpp>

//! Class which stores the grid values of the initial vars in separate arrays.
//! This does not interface with boost and should only be constructed once the
//! integration has finished.

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarSolution
{
public:
    //! Constructor which accepts the arrays outputted by the solution observer
    //! as inputs and puts the data into its own internal arrays.
    BosonStarSolution(initial_data_t<initial_state_t> &a_initial_var_arrays,
    initial_data_t<double> &a_initial_grid);

    //! Returns the maximum radius of the solution
    double get_max_radius() const;

    //! Returns the value of psi at the maximum radius
    double get_psi_at_max_radius() const;

    //! Returns the number of roots in Psi
    int get_num_psi_roots() const;

    //! Calculates the complex oscillation frequency (divided by the scalar
    //! mass) of the solution and stores it in m_frequency
    void calculate_frequency();

    //! Returns the frequency (divided by the scalar mass)
    double get_frequency();

    //! Returns m_initial_grid
    initial_data_t<double>& get_grid();

    //! Returns the m_alpha_array
    initial_data_t<double>& get_alpha();

    //! Returns the m_beta_array
    initial_data_t<double>& get_beta();

    //! Returns the m_psi_array
    initial_data_t<double>& get_psi();

    //! Returns the m_Psi_array
    initial_data_t<double>& get_Psi();

private:
    initial_data_t<double> m_initial_grid; //!< array holding the grid radial coordinates
    initial_data_t<double> m_alpha_array = {}; //!< array to hold alpha grid values
    initial_data_t<double> m_beta_array = {}; //!< array to hold beta grid values
    initial_data_t<double> m_psi_array = {}; //!< array to hold psi grid values
    initial_data_t<double> m_Psi_array = {}; //!< array to hold Psi grid values
    int m_num_grid_values; //!< this stores the number of grid values per variable
    int m_num_psi_roots; //!< this stores the number of roots in psi
    double m_frequency; //!< this stores the complex oscillation frequency/m

    //! Function called by constructor to separate the arrays for each of the
    //! initial variables
    void separateArrays(initial_data_t<initial_state_t> &a_initial_var_arrays);

    //! Function to count the number of roots in psi and store in m_num_psi_roots
    void count_num_psi_roots();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */
