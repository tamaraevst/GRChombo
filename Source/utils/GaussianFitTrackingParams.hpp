/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef GAUSSFITPARAMS_HPP_
#define GAUSSFITPARAMS_HPP_

//! A structure for the input params for the boson star
struct GaussFit_params_t
{
    int num_points; // numer of gridpoints used to create boson star
    int do_star_tracking;
    double offset;
    int field_index;
    double search_width;
};

#endif /* GAUSSFITPARAMS_HPP_ */
