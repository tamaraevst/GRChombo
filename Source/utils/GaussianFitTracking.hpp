/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GAUSSIANFITTRACKING_HPP_
#define GAUSSIANFITTRACKING_HPP_


class GaussianFitTracking
{
  public:
    double m_number = 11.3;

    GaussianFitTracking(double a_number):m_number(a_number){}
    ~GaussianFitTracking() {}
    void Shout()
};

#include "GaussianFitTracking.impl.hpp"

#endif /* GAUSSIANFITTRACKING_HPP_ */
