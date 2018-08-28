/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

//! Class which stores the grid values of the initial vars in separate arrays.
//! This does not interface with boost and should only be constructed once the
//! integration has finished.

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarSolution
{
public:
    //! Constructor 1 that counts the number of roots in psi
    BosonStarSolution(initial_data_t<initial_state_t> &a_initial_var_arrays,
    initial_data_t<double> &a_radii);

    //! Constructor 2 that needs to be passed the number of roots in psi
    BosonStarSolution(initial_data_t<initial_state_t> &a_initial_var_arrays,
    initial_data_t<double> &a_radii, int a_num_psi_roots);

    //! Function called by constructor to separate the arrays for each of the
    //! initial variables
    void separateArrays(initial_data_t<initial_state_t> &a_initial_var_arrays);

    //! Function to return the maximum radius of the solution
    double get_max_radius() const;

    //! Function to retun the value of psi at the maximum get_max_radius
    double get_psi_at_max_radius() const;

    //! Function to count the number of roots in psi if this is not passed to
    //! to the constructor
    void count_num_psi_roots();

    //! Returns the number of roots in Psi
    int get_num_psi_roots() const;

    //! Returns m_radii
    initial_data_t<double>& get_grid();

    //! Returns the m_alpha_array
    initial_data_t<double>& get_alpha();

    //! Returns the m_beta_array
    initial_data_t<double>& get_beta();

    //! Returns the m_psi_array
    initial_data_t<double>& get_psi();

    //! Returns the m_Psi_array
    initial_data_t<double>& get_Psi();

protected:
    initial_data_t<double> m_radii; //!< array holding the grid radial coordinates
    initial_data_t<double> m_alpha_array = {}; //!< array to hold alpha grid values
    initial_data_t<double> m_beta_array = {}; //!< array to hold beta grid values
    initial_data_t<double> m_psi_array = {}; //!< array to hold psi grid values
    initial_data_t<double> m_Psi_array = {}; //!< array to hold Psi grid values
    int m_num_grid_values; //!< this stores the number of grid values per variable
    int m_num_psi_roots; //!< this stores the number of roots in psi
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */
