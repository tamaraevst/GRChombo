/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <algorithm>
#include "StarTracker.hpp"
#include "ChomboParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data
#include "DebuggingTools.hpp"
#include "nr3.h"
// #include "GaussJ.hpp"
#include "FitMRQ.hpp"
#include "FitExample.hpp"

double StarTracker::gaussian(double x, double a, double b, double c)
{
    const double z = (x - b) / c;
    return a * exp(-0.5 * z * z);
}

double StarTracker::find_centre(int a_field_index, int num_star, int direction)
{
    double delta;

    std::vector<double> x_coords(m_points);
    std::vector<double> y_coords(m_points);
    std::vector<double> z_coords(m_points);

    std::vector<double> sigma_vector(m_points);
    std::vector<double> a_vector(3);
    std::vector<double> vals(m_points);
    std::vector<double> vals_f(m_points);

    for (int i = 0; i < m_points; i++)
    {
        delta = m_width * (double(i) / double(m_points - 1) - 0.5);

        if (direction == 0 )
        {
            x_coords[i] = m_star_coords[3 * num_star] + delta;
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 1 )
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1] + delta;
            z_coords[i] = m_star_coords[3 * num_star + 2];
        }

        if (direction == 2 )
        {
            x_coords[i] = m_star_coords[3 * num_star];
            y_coords[i] = m_star_coords[3 * num_star + 1];
            z_coords[i] = m_star_coords[3 * num_star + 2] + delta;
        }

        sigma_vector[i] = 1.0;
    }
 
    m_interpolator->refresh();
    InterpolationQuery query(m_points);
    query.setCoords(0, x_coords.data());
    query.setCoords(1, y_coords.data());
    query.setCoords(2, z_coords.data());
    query.addComp(a_field_index, vals.data(), Derivative::LOCAL,
                  VariableType::evolution);
    m_interpolator->interp(query);

    for (int i = 0; i < m_points; i++)
    {
    	vals_f[i] = 1 - vals[i];	
    }
     
    a_vector[0] = 1 - vals[int(m_points/2)];
    a_vector[1] = m_star_coords[3 * num_star];
    a_vector[2] = 1.;

    if (direction == 0 )
    {
        Fitmrq fitmrq1(x_coords, vals_f, sigma_vector, a_vector, fgauss);
        fitmrq1.fit();
        return fitmrq1.a[1];
    }

    if (direction == 1 )
    {
        Fitmrq fitmrq1(y_coords, vals_f, sigma_vector, a_vector, fgauss);
        fitmrq1.fit();
        return fitmrq1.a[1];
    }

    if (direction == 2 )
    {
        Fitmrq fitmrq1(z_coords, vals_f, sigma_vector, a_vector, fgauss);
        fitmrq1.fit();
        return fitmrq1.a[1];
    }
}

void StarTracker::update_star_centres(int a_field_index)
{
    double starA_0 = find_centre(a_field_index, 0, 0);
    m_star_coords[0] = starA_0;
    double starA_1 = find_centre(a_field_index, 0, 1);
    m_star_coords[1] = starA_1;
    double starA_2 = find_centre(a_field_index, 0, 2);
    m_star_coords[2] = starA_2;

    double starB_0 = find_centre(a_field_index, 1, 0);
    m_star_coords[3] = starB_0;
    double starB_1 = find_centre(a_field_index, 1, 1);
    m_star_coords[4] = starB_1;
    double starB_2 = find_centre(a_field_index, 1, 2);
    m_star_coords[5] = starB_2;
}

void StarTracker::write_to_dat(std::string a_filename, double a_dt,
                               double a_time, double a_restart_time,
                               bool a_first_step)
{
    double eps = 10e-8;
    SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                 SmallDataIO::APPEND, a_first_step);

    if (a_time > a_restart_time + eps)
        star_centre_file.remove_duplicate_time_data();

    std::vector<string> header_line(3. * m_num_stars);

    for (int n = 0; n < m_num_stars; n++)
    {
        header_line[3 * n] = "x" + to_string(n+1);
        header_line[3 * n + 1] = "y" + to_string(n+1);
        header_line[3 * n + 2] = "z" + to_string(n+1);
    }

    if (a_time == 0.)
    {
        star_centre_file.write_header_line(header_line);
    }

    star_centre_file.write_time_data_line(m_star_coords);
}

// void StarTracker::get_star_centres(std::vector<double> &a_centre)
// {
//     int i_max = m_num_stars * 3;
//     a_centre.resize(i_max);
//     for (int i = 0; i < i_max; i++)
//     {
//         a_centre[i] = m_star_coords[i];
//     }
// }

// read a data line from the previous timestep
void StarTracker::read_old_centre_from_dat(std::string a_filename, double a_dt,
                                           double a_time, double a_restart_time,
                                           bool a_first_step)
{
    std::vector<double> data_line;
    if (a_time > a_dt / 3.)
    {
        SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                     SmallDataIO::READ, a_first_step);
        star_centre_file.get_specific_data_line(data_line, a_time - a_dt);

        bool length_match = data_line.size() == m_num_stars * CH_SPACEDIM;

        if (length_match)
        {
            for (int i = 0; i < data_line.size(); i++)
            {
                m_star_coords[i] = data_line[i];
            }
        }
        else
        {
            for (int i = 0; i < m_star_coords.size(); i++)
            {
                m_star_coords[i] = NAN;
            }
            std::cout
                << "Array Size Mismatch While Loading Star Centre From File ! "
                << std::endl;
        }
    }
}

// pass this an empty std::vector and it will resize and fill with centre values
// of field
// void StarTracker::get_field_value_at_centres(
//     int a_field_index, std::vector<double> &a_out_data,
//     AMRInterpolator<Lagrange<4>> *a_interpolator)
// {
//     std::cout << "Getting field centre values" << std::endl;
//     bool gone_NAN = false;
//     a_out_data.resize(m_num_stars);

//     // detect if gone stars centres have gone nan
//     for (int i = 0; i < m_num_stars * CH_SPACEDIM; i++)
//     {
//         if (std::isnan(m_star_coords[i]))
//         {
//             gone_NAN = true;
//         }
//     }

//     if (gone_NAN == false)
//     {

//         std::vector<double> x(m_num_stars);
//         std::vector<double> y(m_num_stars);
//         std::vector<double> z(m_num_stars);
//         std::vector<double> f(m_num_stars);

//         for (int n = 0; n < m_num_stars; n++)
//         {
//             x[n] = m_star_coords[n * CH_SPACEDIM];
//             y[n] = m_star_coords[n * CH_SPACEDIM + 1];
//             z[n] = m_star_coords[n * CH_SPACEDIM + 2];
//         }

//         a_interpolator->refresh();
//         InterpolationQuery query(m_num_stars);
//         query.setCoords(0, x.data());
//         query.setCoords(1, y.data());
//         query.setCoords(2, z.data());
//         query.addComp(a_field_index, f.data());
//         a_interpolator->interp(query);

//         for (int n = 0; n < m_num_stars; n++)
//         {
//             a_out_data[n] = f[n];
//         }
//     }
//     else
//     {
//         for (int n = 0; n < m_num_stars; n++)
//         {
//             a_out_data[n] = NAN;
//         }
//     }
// }
