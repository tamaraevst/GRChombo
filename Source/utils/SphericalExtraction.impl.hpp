/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SPHERICALEXTRACTION_HPP_)
#error "This file should only be included through SphericalExtraction.hpp"
#endif

#ifndef SPHERICALEXTRACTION_IMPL_HPP_
#define SPHERICALEXTRACTION_IMPL_HPP_

//! Set up and execute the interpolation query
inline void SphericalExtraction::execute_query(
    AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    if (a_interpolator == nullptr)
    {
        MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }

    m_interp_var.resize(m_num_points * m_params.num_extraction_radii);
    m_interp_x.resize(m_num_points * m_params.num_extraction_radii);
    m_interp_y.resize(m_num_points * m_params.num_extraction_radii);
    m_interp_z.resize(m_num_points * m_params.num_extraction_radii);

    // Work out the coordinates
    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;
            m_interp_x[iradius * m_num_points + idx] =
                m_params.extraction_center[0] +
                m_params.extraction_radii[iradius] * sin(theta) * cos(phi);
            m_interp_y[iradius * m_num_points + idx] =
                m_params.extraction_center[1] +
                m_params.extraction_radii[iradius] * sin(theta) * sin(phi);
            m_interp_z[iradius * m_num_points + idx] =
                m_params.extraction_center[2] +
                m_params.extraction_radii[iradius] * cos(theta);
        }
    }

    // set up the query
    InterpolationQuery query(m_num_points * m_params.num_extraction_radii);
    query.setCoords(0, m_interp_x.data())
        .setCoords(1, m_interp_y.data())
        .setCoords(2, m_interp_z.data())
        .addComp(m_extraction_comp, m_interp_var.data());

    // submit the query
    a_interpolator->interp(query);
}

//! Write out the result of the extraction in phi and theta at each timestep
inline void SphericalExtraction::write_extraction(
                std::string a_file_prefix) const
{
    CH_TIME("SphericalExtraction::write_extraction");
    SmallDataIO extraction_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW);

    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        // Write headers
        std::vector<std::string> header1_strings = {
            "time = " + std::to_string(m_time) + ",",
            "r = " + std::to_string(m_params.extraction_radii[iradius])};
        extraction_file.write_header_line(header1_strings, "");
        std::vector<std::string> components = {
            UserVariables::variable_names[m_extraction_comp]};
        std::vector<std::string> coords = {"theta", "phi"};
        extraction_file.write_header_line(components, coords);

        // Now the data
        for (int idx = iradius * m_num_points;
             idx < (iradius + 1) * m_num_points; ++idx)
        {
            int itheta =
                (idx - iradius * m_num_points) / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;

            extraction_file.write_data_line({m_interp_var[idx]}, {theta, phi});
        }
        extraction_file.line_break();
    }
}

//! integrate over a spherical shell
inline std::vector<double> SphericalExtraction::integrate_surface() const
{
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    std::vector<double> integral(m_params.num_extraction_radii, 0.);

    // only rank 0 does the integral
    if (rank == 0)
    {
        // integrate the values over the sphere (normalised by r^2) for each
        // radius, assumes spacings constant, uses trapezium rule for phi and
        // rectangles for theta  note we don't have to fudge the end points for
        // phi because the function is periodic  and so the last point (implied
        // but not part of vector) is equal to the first point
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            for (int iphi = 0; iphi < m_params.num_points_phi; ++iphi)
            {
                //double phi = iphi * 2 * M_PI / m_params.num_points_phi;
                double inner_integral = 0.;
                for (int itheta = 0; itheta < m_params.num_points_theta; itheta++)
                {
                    double theta = (itheta + 0.5) * m_dtheta;
                    int idx = iradius * m_num_points +
                              itheta * m_params.num_points_phi + iphi;
                    double f_theta_phi = m_interp_var[idx] * sin(theta);
                    inner_integral += m_dtheta * f_theta_phi;
                }
                integral[iradius] += m_dphi * inner_integral;
            }
        }
    }
    return integral;
}

//! Write out calculated value of integral
inline void SphericalExtraction::write_integral(
            const std::vector<double> a_integral, std::string a_filename) const
{
    CH_TIME("SphericalExtraction::write_integral");
    // open file for writing
    SmallDataIO integral_file(a_filename, m_dt, m_time, m_restart_time,
                              SmallDataIO::APPEND);

    // remove any duplicate data if this is a restart
    // note that this only does something if this is the first timestep after
    // a restart
    integral_file.remove_duplicate_time_data();

    // need to write headers if this is the first timestep
    if (m_time == m_dt)
    {
        // make header strings
        std::vector<std::string> header1_strings(m_params.num_extraction_radii);
        std::vector<std::string> header2_strings(m_params.num_extraction_radii);
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            header1_strings[iradius] = "ADM mass";
            header2_strings[iradius] =
                std::to_string(m_params.extraction_radii[iradius]);
        }

        // write headers
        integral_file.write_header_line(header1_strings);
        integral_file.write_header_line(header2_strings, "r = ");
    }

    // write data
    integral_file.write_time_data_line(a_integral);
}

#endif /* SPHERICALEXTRACTION_IMPL_HPP_ */
