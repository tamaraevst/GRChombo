/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(GAUSSIANFITTRACKING_HPP_)
#error "This file should only be included through GaussianFitTracking.hpp"
#endif

#ifndef GAUSSIANFITTRACKING_IMPL_HPP_
#define GAUSSIANFITTRACKING_IMPL_HPP_

// effectively the main() function
void GaussianFitTracking::do_star_tracking(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    // if t>0 find the previous star centre(s) and decide if theyre usable
    read_old_centre_from_dat();
    check_if_star_positions_are_good();
    get_data(a_interpolator);
    find_centres();
    check_for_BH(a_interpolator);
    write_to_dat();

}

// pass this an empty std::vector and it will resize and fill with centre values of field
void GaussianFitTracking::get_field_value_at_centres(int a_field_index,
                                          std::vector<double> &a_out_data,
                                          AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    a_out_data.resize(1); //make it length 1

    std::vector<double> x(1);
    std::vector<double> y(1);
    std::vector<double> z(1);
    std::vector<double> f(1);

    x[0] = m_centre[0];
    y[0] = m_centre[1];
    z[0] = m_centre[2];

    if ( std::isnan(x[0]) || std::isnan(y[0]) || std::isnan(z[0]))
    {
        a_out_data[0] = NAN;
    }
    else
    {
        a_interpolator->refresh();
        InterpolationQuery query(1);
        query.setCoords(0, x.data());
        query.setCoords(1, y.data());
        query.setCoords(2, z.data());
        query.addComp(a_field_index, f.data());
        a_interpolator->interp(query);
        a_out_data[0] = f[0];
    }

    if (!m_params_GaussFit.track_both_centres) return;

    x[0] = m_centre[3];
    y[0] = m_centre[4];
    z[0] = m_centre[5];

    if ( std::isnan(x[0]) || std::isnan(y[0]) || std::isnan(z[0]))
    {
        a_out_data.push_back(NAN);
    }
    else
    {
        a_interpolator->refresh();
        InterpolationQuery query(1);
        query.setCoords(0, x.data());
        query.setCoords(1, y.data());
        query.setCoords(2, z.data());
        query.addComp(a_field_index, f.data());
        a_interpolator->interp(query);
        a_out_data.push_back(f[0]);
    }
}

/* after the first timestep, decides if the star's positions are good
for use in next tracking */
void GaussianFitTracking::check_if_star_positions_are_good()
{
    if (m_time == 0)
    {
        return;
    }
    // complain if any star has gotten to the edge of the grid
    for (int i=0; i<m_N; i++)
    {
        if (m_old_centre[i]<0.5 || m_old_centre[i] > m_L-0.5)
        {
            return;
        }
    }

    return;
}

// read a data line from the previous timestep
void GaussianFitTracking::read_old_centre_from_dat()
{
    std::vector<double> data_line;
    if (m_time > 10e-8)
    {
        SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                      m_restart_time,
                                      SmallDataIO::READ,
                                      m_first_step);
        star_centre_file.get_specific_data_line(data_line, m_time-m_dt);
    }
    for (int i=0; i<data_line.size(); i++)
    {
        m_old_centre[i] = data_line[i];
    }
}

// loads previously saved data from dat file into m_old_centre arrays
// not called anymore, but left incase its needed for other purposes
void GaussianFitTracking::read_old_centre_from_dat_manually()
{
    if (m_time > 0.)
    {
        int num = 1+m_N; //number of columns of dat file
        std::string messages[num];
        double nums[num];
        std::string message;
        std::ifstream file;
        file.open(m_filename + ".dat", std::ifstream::in);

        file.seekg(-3, file.end); //goes 3 characters before the end of datfile
        //int length = file.tellg(); // character number at end of datfile
        char c;

        // loop marches backwards untill the start of the last line
      while (file.good())
        {
            c = file.get();
            file.unget();
            file.unget();
            if (c=='\n')
            {
              break;
            }
        }
        // go foreward 2 places to miss the '\n' char
        c = file.get();
        c = file.get();

        // walks through the last line splitting the numbers into an array
        bool found_a_number = false; // true if current position in line is not blank or '/n'
        for (int i=0; i<num; i++)
        {
            while (file.good())
            {
                c = file.get();
                message += c;
                if (c!=' ' && c!= '\n' )
                {
                    found_a_number = true;
                }
                if (c==' ' && found_a_number)
                {
                    found_a_number = false;
                    messages[i] = message;
                    nums[i] = std::stod(message);
                    message = " ";
                    break;
                }
                if (c=='\n')
                {
                    messages[i] = message;
                    nums[i] = std::stod(message);
                    break;
                }
            }
        }

        // load the strings read off into the (double) arrays for use
        for (int i=0; i<num; i++)
        {
            if (i!=0)
            {
                m_old_centre[i-1] = nums[i];
            }
        }
        file.close();
    }
}

// uses interpolator to grab field values at desired positions
// also deals with if previous positions were NAN
// also determines the positions
void GaussianFitTracking::get_data(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    double dummy_vals1[3*m_number];
    double dummy_vals2[3*m_number];
    //setup initial star centres from params file
    if (m_time == 0.)
    {
        for (int i=0; i<m_N; i++)
        {
            m_old_centre[i] = 0.5*m_L + m_params_GaussFit.track_centres[i];
        }
    }
    // setup the coordinates required for labeling interpolation sites
    // set of points label 3 mutually orthogonal lines, meeting at centres
    //does this for each star!
    for (int i = 0; i < 3*m_number; i++)
    {
        m_array_x[i] = m_old_centre[0];
        m_array_y[i] = m_old_centre[1];
        m_array_z[i] = m_old_centre[2];
    }
    for (int i = 0; i < m_number; i++)
    {
        m_array_x[i] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_y[i + m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_z[i + 2*m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
    }

    if ( std::isnan(m_old_centre[0]) || std::isnan(m_old_centre[1]) ||
                                                   std::isnan(m_old_centre[2]))
    {
        for (int i=0; i<3*m_number; i++)
        {
            m_vals[i] = NAN;
        }
    }
    else
    {
        // setup interpolator
        a_interpolator->refresh();
        InterpolationQuery query(3*m_number);
        query.setCoords(0, m_array_x);
        query.setCoords(1, m_array_y);
        query.setCoords(2, m_array_z);
        query.addComp(m_field_index, dummy_vals1);
        // submit the query
        a_interpolator->interp(query);

        for (int i=0; i<3*m_number; i++)
        {
            m_vals[i] = dummy_vals1[i];
        }
    }

    if (!m_params_GaussFit.track_both_centres) return;

    for (int i = 0; i < 3*m_number; i++)
    {
        m_array_x[i] = m_old_centre[3];
        m_array_y[i] = m_old_centre[4];
        m_array_z[i] = m_old_centre[5];
    }
    for (int i = 0; i < m_number; i++)
    {
        m_array_x[i] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_y[i + m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_z[i + 2*m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
    }

    if ( std::isnan(m_old_centre[3]) || std::isnan(m_old_centre[4]) ||
                                                   std::isnan(m_old_centre[5]))
    {
        for (int i=3*m_number; i<6*m_number; i++)
        {
            m_vals[i] = NAN;
        }
    }
    else
    {
        // setup interpolator
        a_interpolator->refresh();
        InterpolationQuery query(3*m_number);
        query.setCoords(0, m_array_x);
        query.setCoords(1, m_array_y);
        query.setCoords(2, m_array_z);
        query.addComp(m_field_index, dummy_vals2);
        // submit the query
        a_interpolator->interp(query);

        for (int i=0; i<3*m_number; i++)
        {
            m_vals[i + 3*m_number] = dummy_vals2[i];
        }
    }
}


//calculates the centre of the stars
void GaussianFitTracking::find_centres()
{
    //calculate first star's centre
    if (true)
    {
        double m_integral_x=0., m_integral_y=0., m_integral_z=0.;
        double m_weighted_integral_x=0., m_weighted_integral_y=0.;
        double m_weighted_integral_z=0.;
        for (int i = 0; i < m_number; i++)
        {
            m_integral_x += m_vals[i];
            m_integral_y += m_vals[i+m_number];
            m_integral_z += m_vals[i+2*m_number];
            m_weighted_integral_x += m_vals[i]*m_array_x[i];
            m_weighted_integral_y += m_vals[i+m_number]*m_array_y[i+m_number];
            m_weighted_integral_z += m_vals[i+2*m_number]*m_array_z[i+2*m_number];
        }
        m_centre[0] = m_weighted_integral_x/m_integral_x;
        m_centre[1] = m_weighted_integral_y/m_integral_y;
        m_centre[2] = m_weighted_integral_z/m_integral_z;
    }
    //calculate second star's centre
    if (m_params_GaussFit.track_both_centres)
    {
        double m_integral_x=0., m_integral_y=0., m_integral_z=0.;
        double m_weighted_integral_x=0., m_weighted_integral_y=0.;
        double m_weighted_integral_z=0.;
        for (int i = 0; i < m_number; i++)
        {
            m_integral_x += m_vals[i+3*m_number];
            m_integral_y += m_vals[i+4*m_number];
            m_integral_z += m_vals[i+5*m_number];
            m_weighted_integral_x += m_vals[i+3*m_number]*m_array_x[i];
            m_weighted_integral_y += m_vals[i+4*m_number]*m_array_y[i+m_number];
            m_weighted_integral_z += m_vals[i+5*m_number]*m_array_z[i+2*m_number];
        }
        m_centre[3] = m_weighted_integral_x/m_integral_x;
        m_centre[4] = m_weighted_integral_y/m_integral_y;
        m_centre[5] = m_weighted_integral_z/m_integral_z;
    }
}

void GaussianFitTracking::write_to_dat()
{
    double eps = 10e-8;
    SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                  m_restart_time,
                                  SmallDataIO::APPEND,
                                  m_first_step);
    if (m_time > m_restart_time + eps) star_centre_file.remove_duplicate_time_data();

    if (m_params_GaussFit.track_both_centres)
    {
        if (m_time == 0.)
        {
            star_centre_file.write_header_line({"Star 1 Centre x","Star 1 Centre y",
                                                "Star 1 Centre z","Star 2 Centre x",
                                                "Star 2 Centre y","Star 2 Centre z"});
        }
        star_centre_file.write_time_data_line({m_centre[0],m_centre[1],m_centre[2],
                                               m_centre[3],m_centre[4],m_centre[5]});
    }
    else
    {
        if (m_time == 0.)
        {
            star_centre_file.write_header_line({"Star Centre x","Star Centre y",
                                                "Star Centre z"});
        }
        star_centre_file.write_time_data_line({m_centre[0],m_centre[1],m_centre[2]});
    }
}

void GaussianFitTracking::check_for_BH(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    std::vector<double> central_values; // will be length 1 or 2 depending on num of stars
    get_field_value_at_centres(0,central_values,a_interpolator);
    for (int i=0; i<central_values.size(); i++)
    {
        if (central_values[i]<m_BH_cutoff && !isnan(central_values[i]))
        {
            m_BH_formed.push_back(true);
            for (int j=0; j<3; j++) m_centre[3*i + j] = NAN;
        }
        else
        {
            m_BH_formed.push_back(false);
        }
    }
}

void GaussianFitTracking::get_BH_centres(std::vector<double> &a_out_data)
{
    for (int i=0; i<m_BH_formed.size(); i++)
    {
        for (int j=0; j<3; j++)
        {
            if (m_BH_formed[i])
            {
                a_out_data.push_back(m_old_centre[j+3*i]);
            }
            else
            {
                a_out_data.push_back(NAN);
            }
        }
    }
    /*m_test_file << m_time << ": ";
    for (int i=0; i < m_N; i++)
    {
        m_test_file << a_out_data[i] << ", ";
    }
    m_test_file << std::endl;*/
}

#endif /* GAUSSIANFITTRACKING_IMPL_HPP_ */
