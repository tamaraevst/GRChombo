/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBS_HPP_
#define BINARYBS_HPP_

#include "BosonStar.hpp"

//! Class for initial conditions of two superposed self interacting boson stars
class BinaryBS
{
protected:
    double m_dx;
    const BosonStar m_BosonStar1;
    const BosonStar m_BosonStar2;
    bool m_identical; // if true, then both stars have the same spherically
                    // symmetric profile
    int m_verbosity;

public:
    //! The constructor
    BinaryBS(BosonStar_params_t a_params_BosonStar1,
             BosonStar_params_t a_params_BosonStar2,
             Potential::params_t a_params_potential,
             double a_G_Newton, double a_dx, int a_verbosity);

    //! Computes the spherically symmetric profile for each boson star
    void compute_profiles();

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;
};

#include "BinaryBS.impl.hpp"

#endif /* BINARYBS_HPP_ */
