/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDSCALARS_HPP_
#define MODIFIEDSCALARS_HPP_

#include "simd.hpp"

//This is for switches of any of the non-GR scalars we wish to introduce.

class ModifiedScalars
{
  public:
    struct params_t
    {
         bool csswitch;
         bool gbswitch;
    };

  protected:
    const params_t m_params;

  public:
    //! The constructor
    ModifiedScalars(params_t a_params) : m_params(a_params) {}


};

#endif /* MODIFIEDSCALARS_HPP_ */
