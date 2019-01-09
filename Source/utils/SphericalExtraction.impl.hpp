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
    AMRInterpolator<Lagrange<4>> *m_interpolator) const
{
    if (m_interpolator == nullptr)
    {
        MayDay::Error("Interpolator has not been initialised in GRAMR class.");
    }

    // Work out the coordinates
    for (int idx = 0; idx < m_num_points; ++idx)
    {
        int itheta = idx / m_params.num_points_phi;
        int iphi = idx % m_params.num_points_phi;
        // don't put a point at z = 0
        double theta = (itheta + 0.5) * m_dtheta;
        double phi = iphi * m_dphi;
        m_interp_x[idx] = m_params.extraction_center[0] +
                          m_params.extraction_radius * sin(theta) * cos(phi);
        m_interp_y[idx] = m_params.extraction_center[1] +
                          m_params.extraction_radius * sin(theta) * sin(phi);
        m_interp_z[idx] = m_params.extraction_center[2] +
                          m_params.extraction_radius * cos(theta);
    }

    // set up the query
    InterpolationQuery query(m_num_points);
    query.setCoords(0, m_interp_x.get())
        .setCoords(1, m_interp_y.get())
        .setCoords(2, m_interp_z.get())
        .addComp(m_extraction_comp, m_state_ptr.get());

    // submit the query
    m_interpolator->interp(query);
}

//! Write out the result of the extraction in phi and theta at each timestep
inline void SphericalExtraction::write_extraction(string file_prefix) const
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
        int step = m_time / m_dt;
        string file_str = file_prefix + std::to_string(step);
        char comp_str[20];
        sprintf(comp_str, UserVariables::variable_names[m_extraction_comp]);

        // write out complete data to a separate file at each step
        std::ofstream outfile;
        outfile.open(file_str);
        if (!outfile.is_open())
        {
            MayDay::Error(
                "in SphericalExtraction::error opening output file for "
                "extraction points");
        }

        // header data
        outfile << "# time : " << m_time << endl;
        outfile << "#" << std::setw(19) << "theta";
        outfile << std::setw(20) << "phi";
        outfile << std::setw(20) << comp_str << endl;

        // Now the data
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            int itheta = idx / m_params.num_points_phi;
            int iphi = idx % m_params.num_points_phi;
            // don't put a point at z = 0
            double theta = (itheta + 0.5) * m_dtheta;
            double phi = iphi * m_dphi;
            outfile << std::setw(20) << theta;
            outfile << std::setw(20) << phi;
            outfile << std::setw(20) << std::setprecision(9) << m_state_ptr[idx]
                    << endl;
        }
        outfile.close();
    }
}

//! integrate over a spherical shell
inline double SphericalExtraction::integrate_surface() const
{
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
#else
    rank = 0;
#endif
    std::vector<double> integrand;
    double integral = 0.;

    // only rank 0 does the integral
    if (rank == 0)
    {
        for (int idx = 0; idx < m_num_points; ++idx)
        {
            // setup the integrand for next stage
            double x = m_interp_x[idx] - m_params.extraction_center[0];
            double y = m_interp_y[idx] - m_params.extraction_center[1];
            double z = m_interp_z[idx] - m_params.extraction_center[2];
            integrand[idx] = m_state_ptr[idx];
        }
        // integrate the values over the sphere (normalised by r^2)
        // assumes spacings constant, uses trapezium rule for phi and rectangles for
        // theta  note we don't have to fudge the end points for phi because the
        // function is periodic  and so the last point (implied but not part of
        // vector) is equal to the first point
        for (int iphi = 0; iphi < m_params.num_points_phi; ++iphi)
        {
            double phi = iphi * 2 * M_PI / m_params.num_points_phi;
            double inner_integral = 0;
            for (int itheta = 0; itheta < m_params.num_points_theta; itheta++)
            {
                double theta = (itheta + 0.5) * m_dtheta;
                int idx = itheta * m_params.num_points_phi + iphi;
                double f_theta_phi = integrand[idx] * sin(theta);
                inner_integral += m_dtheta * f_theta_phi;
            }
            integral += m_dphi * inner_integral;
        }
    }
    return integral;
}

//! Write out calculated value of integral
inline void SphericalExtraction::write_integral(double integral,
                                                string file_name) const
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
        // and append the integral output to a file - one file for whole run
        int step = m_time / m_dt;
        string file_str = file_name;
        std::ofstream outfile2;
        if (m_time == m_dt)
        {
            outfile2.open(file_str);
        }
        else
        {
            outfile2.open(file_str, std::ios_base::app);
        }
        if (!outfile2.is_open())
        {
            MayDay::Error(
                "in SphericalExtraction::error opening output file for "
                "extraction integral");
        }

        // Header data at first timestep
        if (m_time == m_dt)
        {
            outfile2 << "#" << std::setw(19) << "m_time";
            outfile2 << std::setw(20) << "integral" << std::endl;
        }
        outfile2 << std::setw(20) << m_time;
        outfile2 << std::setw(20) << integral << std::setprecision(9)
                 << std::endl;
        outfile2.close();
    }
}

#endif /* SPHERICALEXTRACTION_IMPL_HPP_ */
