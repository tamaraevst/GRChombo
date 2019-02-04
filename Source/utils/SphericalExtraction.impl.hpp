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
    AMRInterpolator<Lagrange<4>> *a_interpolator) const
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
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif

    // only rank 0 does the write out
    if (rank == 0)
    {
        // set up file names and component names
        int step = std::round(m_time / m_dt);
        std::string file_str = a_file_prefix + std::to_string(step);
        std::string comp_str = UserVariables::variable_names[m_extraction_comp];

        // write out complete data to a separate file at each step
        std::ofstream outfile;
        outfile.open(file_str);
        if (!outfile)
        {
            MayDay::Error(
                "SphericalExtraction::write_extraction: error opening output "
                "file");
        }
        // header data
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            outfile << "# time : " << m_time << ", r = ";
            outfile << m_params.extraction_radii[iradius] << "\n";
            outfile << "#" << std::setw(11) << "theta";
            outfile << std::setw(12) << "phi";
            outfile << std::setw(20) << comp_str << "\n";

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
                outfile << std::fixed << std::setprecision(7);
                outfile << std::setw(12) << theta;
                outfile << std::setw(12) << phi;
                outfile << std::scientific << std::setprecision(10);
                outfile << std::setw(20) << m_interp_var[idx] << "\n";
            }
            outfile << "\n\n";
        }
        outfile.close();
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
inline void SphericalExtraction::write_integral(const std::vector a_integral,
                                                std::string a_filename) const
{
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    // only rank 0 does the write out
    if (rank == 0)
    {
        // overwrite file if this is the first timestep, otherwise append.
        std::ofstream outfile;
        if (m_time == m_dt)
        {
            outfile.open(a_filename);
        }
        else
        {
            outfile.open(a_filename, std::ios_base::app);
        }
        if (!outfile)
        {
            MayDay::Error(
                "SphericalExtraction::write_integral: error opening output "
                "file");
        }

        if (m_time == m_dt)
        {
            outfile << "#" << std::setw(9) << "time";
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                outfile << std::setw(20) << "integral";
            }
            outfile << "\n";
            outfile << "#" << std::setw(9) << "r =";
            for (int iradius = 0;
                 iradius < m_params.num_extraction_radii; ++iradius)
            {
                outfile << std::setw(20) << m_params.extraction_radii[iradius];
            }
            outfile << "\n";
        }


        outfile << std::fixed << std::setw(10) << m_time;
        outfile << std::scientific << std::setprecision(10);
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            outfile << std::setw(20) << a_integral[iradius];
        }
        outfile << std::endl;
        outfile.close();
    }
}

#endif /* SPHERICALEXTRACTION_IMPL_HPP_ */
