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
    double central_amplitude_CSF; //!< Central amplitude of the star
    double phase;
    int eigen; // radial eigenstate of the boson star (0=ground)
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
};

#endif /* BOSONSTARPARAMS_HPP_ */
