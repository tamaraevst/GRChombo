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
    int field_index;
    int AMR_level;
    double search_width;
    bool track_both_centres;
    double track_min_separation;
    double BH_cutoff;
    std::array<double, CH_SPACEDIM> track_centre; //!< coordinates of the centres of the stars
    std::array<double, 2*CH_SPACEDIM> track_centres; //!< coordinates of the centres of the stars
};

#endif /* GAUSSFITPARAMS_HPP_ */