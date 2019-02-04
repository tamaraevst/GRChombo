/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MASSEXTRACTION_HPP_
#define MASSEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
#include "UserVariables.hpp" // Needs c_Madm

class MassExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    MassExtraction(extraction_params_t a_params, double a_dt, double a_time)
        : SphericalExtraction(c_Madm, a_params, a_dt, a_time) {}

    //! Extract the mass
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        CH_TIME("MassExtraction::execute_query");
        SphericalExtraction::execute_query(a_interpolator);
        std::vector<double> Madm = SphericalExtraction::integrate_surface();
        SphericalExtraction::write_integral(Madm, "ADMmass.dat");
    }
};

#endif /* MASSEXTRACTION_HPP_ */
