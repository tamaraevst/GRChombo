/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(PUNCTURETRACKER_HPP_)
#error "This file should only be included through PunctureTracker.hpp"
#endif

#ifndef PUNCTURETRACKER_IMPL_HPP_
#define PUNCTURETRACKER_IMPL_HPP_

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures(
    GRAMR &a_gramr,
    std::vector<std::array<double, CH_SPACEDIM>> initial_puncture_coords) const
{
    // first set the puncture data, initial shift is always zero
    std::vector<std::array<double, CH_SPACEDIM>> initial_shift;
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR1(i) { initial_shift[ipuncture][i] = 0.0; }
    }
    a_gramr.set_puncture_data(initial_puncture_coords, initial_shift);

    // now the write out
    SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                               m_restart_time, SmallDataIO::APPEND);
    std::vector<std::string> header1_strings(CH_SPACEDIM * m_num_punctures);
    header1_strings[0] = "time";
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        header1_strings[CH_SPACEDIM * ipuncture + 1] = "x";
        header1_strings[CH_SPACEDIM * ipuncture + 2] = "y";
        header1_strings[CH_SPACEDIM * ipuncture + 3] = "z";
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    std::vector<double> puncture_vector =
        get_puncture_vector(initial_puncture_coords);
    punctures_file.write_time_data_line(puncture_vector);
}

//! Set punctures post restart
void PunctureTracker::read_in_punctures(GRAMR &a_gramr) const
{
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords;
    puncture_coords.resize(m_num_punctures);

    // read them in from the Punctures file at current time m_time
    SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                               m_restart_time, SmallDataIO::READ);
    // NB need to give the get function an empty vector to fill
    std::vector<double> puncture_vector;
    punctures_file.get_specific_data_line(puncture_vector, m_time);
    CH_assert(puncture_vector.size() == m_num_punctures * CH_SPACEDIM);
    punctures_file.remove_duplicate_time_data();

    // convert vector to list of coords
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_coords[ipuncture] = {
            puncture_vector[ipuncture * CH_SPACEDIM + 0],
            puncture_vector[ipuncture * CH_SPACEDIM + 1],
            puncture_vector[ipuncture * CH_SPACEDIM + 2]};
    }

    // set the coordinates and get the current shift
    std::vector<std::array<double, CH_SPACEDIM>> current_shift;
    current_shift.resize(m_num_punctures);
    current_shift = get_interp_shift(a_gramr, puncture_coords);
    a_gramr.set_puncture_data(puncture_coords, current_shift);

    // print out values into pout files
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        pout() << "Puncture " << ipuncture
               << " restarted at : " << puncture_coords[ipuncture][0] << " "
               << puncture_coords[ipuncture][1] << " "
               << puncture_coords[ipuncture][2] << endl;
        pout() << "at time = " << m_time << endl;
    }
}

//! Execute the tracking and write out
void PunctureTracker::execute_tracking(GRAMR &a_gramr,
                                       const bool write_punctures) const
{
    // get puncture coordinates and old shift value
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords =
        a_gramr.get_puncture_coords();
    std::vector<std::array<double, CH_SPACEDIM>> old_shift =
        a_gramr.get_puncture_shift();
    CH_assert(puncture_coords.size() == m_num_punctures); // sanity check

    // new shift value
    std::vector<std::array<double, CH_SPACEDIM>> new_shift;
    new_shift.resize(m_num_punctures);
    new_shift = get_interp_shift(a_gramr, puncture_coords);

    // update puncture locations using second order update
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR1(i)
        {
            puncture_coords[ipuncture][i] +=
                -0.5 * m_dt *
                (new_shift[ipuncture][i] + old_shift[ipuncture][i]);
        }
    }
    a_gramr.set_puncture_data(puncture_coords, new_shift);

    // print them out
    if (write_punctures)
    {
        SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                                   m_restart_time, SmallDataIO::APPEND);

        // use a vector for the write out
        std::vector<double> puncture_vector =
            get_puncture_vector(puncture_coords);
        punctures_file.write_time_data_line(puncture_vector);
    }
}

//! Use the interpolator to get the value of the shift at
//! given coords
std::vector<std::array<double, CH_SPACEDIM>> PunctureTracker::get_interp_shift(
    GRAMR &a_gramr,
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const
{
    // interpolated shift value
    std::vector<std::array<double, CH_SPACEDIM>> interp_shift;
    interp_shift.resize(m_num_punctures);

    // refresh interpolator
    a_gramr.m_interpolator->refresh();

    // set up shift and coordinate holders
    std::vector<double> interp_shift1(m_num_punctures);
    std::vector<double> interp_shift2(m_num_punctures);
    std::vector<double> interp_shift3(m_num_punctures);
    std::vector<double> interp_x(m_num_punctures);
    std::vector<double> interp_y(m_num_punctures);
    std::vector<double> interp_z(m_num_punctures);

    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_x[ipuncture] = puncture_coords[ipuncture][0];
        interp_y[ipuncture] = puncture_coords[ipuncture][1];
        interp_z[ipuncture] = puncture_coords[ipuncture][2];
    }

    // setup query
    InterpolationQuery query(m_num_punctures);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data())
        .addComp(c_shift3, interp_shift3.data());

    // engage!
    a_gramr.m_interpolator->interp(query);

    // put the shift values into an array
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_shift[ipuncture] = {interp_shift1[ipuncture],
                                   interp_shift2[ipuncture],
                                   interp_shift3[ipuncture]};
    }
    return interp_shift;
}

//! get a vector of the puncture coords - used for write out
std::vector<double> PunctureTracker::get_puncture_vector(
    std::vector<std::array<double, CH_SPACEDIM>> puncture_coords) const
{
    std::vector<double> puncture_vector;
    puncture_vector.resize(m_num_punctures * CH_SPACEDIM);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_vector[ipuncture * CH_SPACEDIM + 0] =
            puncture_coords[ipuncture][0];
        puncture_vector[ipuncture * CH_SPACEDIM + 1] =
            puncture_coords[ipuncture][1];
        puncture_vector[ipuncture * CH_SPACEDIM + 2] =
            puncture_coords[ipuncture][2];
    }
    return puncture_vector;
}

#endif /* PUNCTURETRACKER_IMPL_HPP_ */
