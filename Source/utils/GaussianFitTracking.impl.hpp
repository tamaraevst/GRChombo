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
    m_test_file << "m_time " << m_time << ", m_restart_time " << m_restart_time << std::endl;
    read_old_centre_from_dat();

    star_positions_are_good=true;//calculate_star_position_bool();

    if (star_positions_are_good)
    {
        if (m_params_GaussFit.track_both_centres)
        {
            get_data_2(a_interpolator);
            find_centres();
            write_to_dat_2();
        }
        else
        {
            get_data(a_interpolator);
            find_centres();
            write_to_dat();
        }
    }
    read_old_centre_from_dat_2();
}

/* after the first timestep, decides if the star's positions are good
for use in next tracking */
bool GaussianFitTracking::calculate_star_position_bool()
{
    if (m_time == 0)
    {
        return true;
    }
    // complain if any star has gotten to the edge of the grid
    for (int i=0; i<m_N; i++)
    {
        if (m_old_centre[i]<0.5 || m_old_centre[i] > m_L-0.5)
        {
            return false;
        }
    }
    // complain if two tracked centres get too close to each other
    if (m_params_GaussFit.track_both_centres)
    {
        double star_sep = (m_old_centre[0]-m_old_centre[3])*(m_old_centre[0]-m_old_centre[3])
                   + (m_old_centre[1]-m_old_centre[4])*(m_old_centre[1]-m_old_centre[4])
                   + (m_old_centre[2]-m_old_centre[5])*(m_old_centre[2]-m_old_centre[5]);
        if (star_sep < m_min_separation*m_min_separation)
        {
            return false;
        }
    }
    return true;
}


void GaussianFitTracking::read_old_centre_from_dat_2()
{
    if (!m_first_step)
    {
        std::vector<double> *data_line;
        SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                      m_restart_time,
                                      SmallDataIO::READ,
                                      m_first_step);
        star_centre_file.get_specific_data_line(data_line, m_time-m_dt);
    }

    m_test_file << data_line[0] << ", " << data_line[1] << std::endl;
}

// loads previously saved data from dat file into m_old_centre arrays
void GaussianFitTracking::read_old_centre_from_dat()
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
                //m_test_file << m_time << " " << c << std::endl;
                message += c;
                if (c!=' ' && c!= '\n' )
                {
                    found_a_number = true;
                }
                if (c==' ' && found_a_number)
                {
                    found_a_number = false;
                    messages[i] = message;
                    //m_test_file << m_time << " " << message << std::endl;
                    nums[i] = std::stod(message);
                    message = " ";
                    break;
                }
                if (c=='\n')
                {
                    messages[i] = message;
                    //m_test_file << m_time << " " << message << std::endl;
                    nums[i] = std::stod(message);
                    break;
                }
            }
        }

        // load the strings read off into the (double) arrays for use
        for (int i=0; i<num; i++)
        {
            //m_test_file << "i " << i << " : msgs " << messages[i] << std::endl;
        }
        for (int i=0; i<num; i++)
        {
            //nums[i] = std::stod(messages[i]);
            if (i!=0)
            {
                m_old_centre[i-1] = nums[i];
                //m_test_file << "i " << i << " : m_level " << m_level << ", i=" << i << std::endl;

            }
        }
        file.close();
    }
}

// THIS FUNCTION IS USED FOR TRACKING ONE STAR ONLY
// uses interpolator to grab field values at desired positions
// also determines the positions
void GaussianFitTracking::get_data(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    //setup initial star centres from params file
    if (m_time == 0.)
    {
        for (int i=0; i<m_N; i++)
        {
            m_old_centre[i] = 0.5*m_L + m_params_GaussFit.track_centre[i];
        }
    }
    // setup the coordinates required for labeling interpolation sites
    // set of points label 3 mutually orthogonal lines, meeting at centres
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


    m_test_file << m_old_centre[0] << ", " << m_old_centre[1] << ", "
                                           << m_old_centre[2] << std::endl;
    // setup the interpolator
    a_interpolator->refresh();
    InterpolationQuery query(3*m_number);
    query.setCoords(0, m_array_x);
    query.setCoords(1, m_array_y);
    query.setCoords(2, m_array_z);
    query.addComp(m_field_index, m_vals);
    // submit the query
    a_interpolator->interp(query);
}

// THIS FUNCTION IS USED FOR TRACKING TWO STARS ONLY
// uses interpolator to grab field values at desired positions
// also determines the positions
void GaussianFitTracking::get_data_2(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
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
    for (int i = 3*m_number; i < 6*m_number; i++)
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
        m_array_x[i + 3*m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_y[i + 4*m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
        m_array_z[i + 5*m_number] += 0.5*m_delta*(2*i + 1 - m_number)/m_number;
    }

    m_test_file << m_old_centre[0] << ", " << m_old_centre[1] << ", "
                                           << m_old_centre[2] << std::endl;
    // setup interpolator
    a_interpolator->refresh();
    InterpolationQuery query(m_N*m_number);
    query.setCoords(0, m_array_x);
    query.setCoords(1, m_array_y);
    query.setCoords(2, m_array_z);
    query.addComp(m_field_index, m_vals);
    // submit the query
    a_interpolator->interp(query);
}


//calculates the centre of the stars
void GaussianFitTracking::find_centres()
{
    //calculate first star's centre
    if (true)
    {
        double m_integral_x=0., m_integral_y=0., m_integral_z=0.;
        double m_weighted_integral_x=0., m_weighted_integral_y=0., m_weighted_integral_z=0.;
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
        double m_weighted_integral_x=0., m_weighted_integral_y=0., m_weighted_integral_z=0.;
        for (int i = 0; i < m_number; i++)
        {
            m_integral_x += m_vals[i+3*m_number];
            m_integral_y += m_vals[i+4*m_number];
            m_integral_z += m_vals[i+5*m_number];
            m_weighted_integral_x += m_vals[i+3*m_number]*m_array_x[i+3*m_number];
            m_weighted_integral_y += m_vals[i+4*m_number]*m_array_y[i+4*m_number];
            m_weighted_integral_z += m_vals[i+5*m_number]*m_array_z[i+5*m_number];
        }
        m_centre[3] = m_weighted_integral_x/m_integral_x;
        m_centre[4] = m_weighted_integral_y/m_integral_y;
        m_centre[5] = m_weighted_integral_z/m_integral_z;
    }
}



void GaussianFitTracking::write_to_dat()
{
    SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                  m_restart_time,
                                  SmallDataIO::APPEND,
                                  m_first_step);
    if (m_time > m_restart_time + eps){star_centre_file.remove_duplicate_time_data();}
    if (m_time == 0.)
    {
        star_centre_file.write_header_line({"Star Centre x","Star Centre y","Star Centre z"});
    }
    star_centre_file.write_time_data_line({m_centre[0],m_centre[1],m_centre[2]});
}

void GaussianFitTracking::write_to_dat_2()
{
    SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                  m_restart_time,
                                  SmallDataIO::APPEND,
                                  m_first_step);
    if (m_time > m_restart_time + eps){star_centre_file.remove_duplicate_time_data();}
    if (m_time == 0.)
    {
        star_centre_file.write_header_line({"Star 1 Centre x","Star 1 Centre y",
                                            "Star 1 Centre z","Star 2 Centre x",
                                            "Star 2 Centre y","Star 2 Centre z"});
    }
    star_centre_file.write_time_data_line({m_centre[0],m_centre[1],m_centre[2],
                                           m_centre[3],m_centre[4],m_centre[5]});
}

#endif /* GAUSSIANFITTRACKING_IMPL_HPP_ */
