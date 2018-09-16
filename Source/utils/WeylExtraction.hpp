/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef WEYLEXTRACTION_HPP_
#define WEYLEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "SphericalHarmonics.hpp"
#include "UserVariables.hpp" // Needs c_Weyl_Re etc
#include <fstream>
#include <iostream>

//!  The class allows extraction of the values of the Weyl scalar components on
//!  a spherical shell, and integration over that shell
/*!
   The class allows the user to extract data from the grid for thr Weyl
   components over a spherical shell. The values may then be written to an
   output file, or integrated across the surface.
*/
class WeylExtraction
{
  private:
    //! Params for extraction
    const extraction_params_t m_params;
    const int m_re_comp = c_Weyl4_Re;
    const int m_im_comp = c_Weyl4_Im;
    const double m_dt;
    const double m_time;
    const int m_num_points;
    const double m_dphi;
    const double m_dtheta;

  public:
    //! The constructor
    WeylExtraction(extraction_params_t a_params, double a_dt, double a_time)
        : m_params(a_params), m_dt(a_dt), m_time(a_time),
          m_num_points(m_params.num_points_phi * m_params.num_points_theta),
          m_dphi(2.0 * M_PI / m_params.num_points_phi),
          m_dtheta(M_PI / m_params.num_points_theta)
    {
    }

    //! Destructor
    ~WeylExtraction() {}

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *m_interpolator) const;

  private:
    //! integrate over a spherical shell with given harmonics
    std::array<double, 2> integrate_surface(int es, int el, int em,
                                            const double *m_state_ptr_re,
                                            const double *m_state_ptr_im) const;

    //! Write out calculated value of integral
    void write_integral(std::array<double, 2> integral,
                        const char *filename) const;

    //! Write out the result of the extraction in phi and theta at each timestep
    void write_extraction(char *file_prefix, const double *m_state_ptr_re,
                          const double *m_state_ptr_im) const;
};

#include "WeylExtraction.impl.hpp"

#endif /* WEYLEXTRACTION_HPP_ */
