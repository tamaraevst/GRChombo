/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SPHERICALEXTRACTION_HPP_
#define SPHERICALEXTRACTION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "UserVariables.hpp"
#include <fstream>
#include <iostream>
#include <memory>

//!  The class allows extraction of the values of a single variable at points on
//!  a spherical shell
/*!
     The class allows the user to extract data from the grid for a single
   component over a spherical shell. The values may then be written to an output
   file, or integrated across the surface.
*/
class SphericalExtraction
{
  private:
    //! Params for extraction
    const extraction_params_t m_params;
    const int m_extraction_comp;
    const double m_dt;
    const double m_time;
    const int m_num_points;
    std::unique_ptr<double[]> m_state_ptr{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_x{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_y{new double[m_num_points]};
    std::unique_ptr<double[]> m_interp_z{new double[m_num_points]};
    const double m_dphi;
    const double m_dtheta;

  public:
    //! The constructor
    SphericalExtraction(int a_extraction_comp, extraction_params_t a_params,
                        double a_dt, double a_time)
        : m_extraction_comp(a_extraction_comp), m_params(a_params), m_dt(a_dt),
          m_time(a_time),
          m_num_points(m_params.num_points_phi * m_params.num_points_theta),
          m_dphi(2.0 * M_PI / m_params.num_points_phi),
          m_dtheta(M_PI / m_params.num_points_theta)
    {}

    //! Execute the interpolation query
    void execute_query(AMRInterpolator<Lagrange<4>> *m_interpolator) const;

    //! Write out the result of the extraction in phi and theta at each timestep
    void write_extraction(string file_prefix = "Interpolation") const;

    //! integrate extracted values over a spherical shell
    double integrate_surface() const;

    //! Write out calculated value of intergal
    void write_integral(double integral, string file_prefix = "Integral") const;
};

#include "SphericalExtraction.impl.hpp"

#endif /* SPHERICALEXTRACTION_HPP_ */
