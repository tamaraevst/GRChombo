/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(GAUSSIANFITTRACKING_HPP_)
#error "This file should only be included through GaussianFitTracking.hpp"
#endif

#ifndef GAUSSIANFITTRACKING_IMPL_HPP_
#define GAUSSIANFITTRACKING_IMPL_HPP_

void GaussianFitTracking::do_star_tracking(AMRInterpolator<Lagrange<4>> *a_interpolator)
{

    


    read_old_centre_from_dat();
    get_data(a_interpolator);
    find_centre();
    write_to_dat();
}

void GaussianFitTracking::read_old_centre_from_dat()
{
    if (!m_first_step)
    {
        int num = 4;
        std::string messages[num];
        double nums[num];
        std::string message;
        std::ifstream file;
        file.open(m_filename + ".dat", std::ifstream::in);

        file.seekg(-3, file.end); //goes 3 characters before the end of the document
        int length = file.tellg(); // unused line
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
                    message = " ";
                    break;
                }
                if (c=='\n')
                {
                    messages[i] = message;
                    break;
                }
            }
        }
        for (int i=0; i<num; i++)
        {
            nums[i] = stod(messages[i]);
            if (i!=0)
            {
                m_old_centre[i-1] = nums[i];
                m_test_file << "Time " << m_time << " : msg " << m_old_centre[i-1] << ", i=" << i << std::endl;
            }
        }
        file.close();
    }
}

void GaussianFitTracking::get_data(AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    //grab the interpolated datapoints
    if (m_first_step)
    {
      m_old_centre[0] = 0.5*m_L + m_params_GaussFit.offset;
      m_old_centre[1] = 0.5*m_L;
      m_old_centre[2] = 0.5*m_L;
    }
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

    a_interpolator->refresh();
    InterpolationQuery query(3*m_number);
    query.setCoords(0, m_array_x);
    query.setCoords(1, m_array_y);
    query.setCoords(2, m_array_z);
    query.addComp(m_field_index, m_vals);
    // submit the query
    a_interpolator->interp(query);
}

void GaussianFitTracking::find_centre()
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


void GaussianFitTracking::write_to_dat()
{
    SmallDataIO star_centre_file(m_filename, m_dt, m_time,
                                  m_restart_time,
                                  SmallDataIO::APPEND,
                                  m_first_step);
    star_centre_file.remove_duplicate_time_data();
    if (m_time == 0.)
    {
        star_centre_file.write_header_line({"Star Centre x","Star Centre y","Star Centre z"});
    }
    star_centre_file.write_time_data_line({m_centre[0],m_centre[1],m_centre[2]});
}


#endif /* GAUSSIANFITTRACKING_IMPL_HPP_ */
