/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

//! A structure for the input params for the boson star
struct BosonStar_params_t
{
    int gridpoints; // numer of gridpoints used to create boson star
    double central_amplitude_CSF; //!< Central amplitude of the star
    double phase;
    int eigen; // radial eigenstate of the boson star (0=ground)
    double BS_separation;
    double BS_impact_parameter;
    bool BS_binary;
    bool BS_BH_binary;
    double Newtons_constant;
    double BlackHoleMass;
    double BS_rapidity;
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star

    int gridpoints2;
    int eigen2;
    double Newtons_constant2;
    double central_amplitude_CSF2; //!< Central amplitude of the star
    double BS_rapidity2;
};

#endif /* BOSONSTARPARAMS_HPP_ */
