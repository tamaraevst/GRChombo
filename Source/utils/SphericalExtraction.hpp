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
#include "SmallDataIO.hpp"
#include "UserVariables.hpp"

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
    const double m_restart_time;
    const bool m_called_in_do_analysis;
    const int m_num_points; // number of points per extraction radius
    std::vector<double> m_interp_var;
    std::vector<double> m_interp_x;
    std::vector<double> m_interp_y;
    std::vector<double> m_interp_z;
    const double m_dphi;
    const double m_dtheta;

  public:
    //! The constructor
    SphericalExtraction(int a_extraction_comp, extraction_params_t a_params,
                        double a_dt, double a_time, double a_restart_time = 0.0,
                        bool a_called_in_do_analysis = false)
        : m_params(a_params), m_extraction_comp(a_extraction_comp), m_dt(a_dt),
          m_time(a_time), m_restart_time(a_restart_time),
          m_called_in_do_analysis(a_called_in_do_analysis),
          m_num_points(m_params.num_points_phi * m_params.num_points_theta),
          m_dphi(2.0 * M_PI / m_params.num_points_phi),
          m_dtheta(M_PI / m_params.num_points_theta)
    {}

    //! Execute the interpolation query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator);

    //! Write out the result of the extraction in phi and theta at each timestep
    void write_extraction(std::string a_file_prefix = "Interpolation") const;

    //! integrate extracted values over a spherical shell
    std::vector<double> integrate_surface() const;

    //! Write out calculated value of intergal
    void write_integral(std::vector<double> a_integral,
                        std::string a_filename = "Integral") const;
};

#include "SphericalExtraction.impl.hpp"

#endif /* SPHERICALEXTRACTION_HPP_ */
