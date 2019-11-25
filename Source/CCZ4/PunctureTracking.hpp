/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PUNCTURETRACKING_HPP_
#define PUNCTURETRACKING_HPP_

#include "AMRInterpolator.hpp"
#include "GRAMR.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SimulationParametersBase.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // Needs c_shift etc
//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class PunctureTracker
{
  private:
    //! Params for puncture tracking
    const int m_num_punctures;
    const double m_time;
    const double m_restart_time;
    const double m_dt;
    const std::string m_punctures_filename;

  public:
    //! The constructor
    PunctureTracker(const double a_time, const double a_restart_time,
                    const double a_dt, const std::string a_checkpoint_prefix,
                    const int a_num_punctures = 2)
        : m_num_punctures(a_num_punctures), m_time(a_time),
          m_punctures_filename(a_checkpoint_prefix + "Punctures.txt"),
          m_restart_time(a_restart_time), m_dt(a_dt)
    {
    }

    //! Set punctures post restart
    void read_in_punctures(GRAMR &a_gramr) const
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

        // convert vector to list of coords
        for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
        {
            puncture_coords[ipuncture] = {
                puncture_vector[ipuncture * CH_SPACEDIM + 0],
                puncture_vector[ipuncture * CH_SPACEDIM + 1],
                puncture_vector[ipuncture * CH_SPACEDIM + 2]};
        }

        // set the coordinates
        a_gramr.set_puncture_coords(puncture_coords);

        // print out values into pout files
        pout() << "Punctures restarted at " << puncture_coords[0][0] << " "
               << puncture_coords[0][1] << " " << puncture_coords[0][2] << endl;
        pout() << "and                    " << puncture_coords[1][0] << " "
               << puncture_coords[1][1] << " " << puncture_coords[1][2] << endl;
        pout() << "at time = " << m_time << endl;
    }

    //! Execute the tracking and write out
    //! If optional arguments unset will write on every timestep
    //! on which the tracking is performed
    void execute_tracking(GRAMR &a_gramr, const bool write_punctures = true,
                          const double a_coarse_dt = 0.0) const
    {
        // refresh interpolator
        a_gramr.m_interpolator->refresh();

        // set up shift and coordinate holders
        std::vector<double> interp_shift1(m_num_punctures);
        std::vector<double> interp_shift2(m_num_punctures);
        std::vector<double> interp_shift3(m_num_punctures);
        std::vector<double> interp_x(m_num_punctures);
        std::vector<double> interp_y(m_num_punctures);
        std::vector<double> interp_z(m_num_punctures);

        // assign coordinates
        std::vector<std::array<double, CH_SPACEDIM>> new_puncture_coords =
            a_gramr.get_puncture_coords();
        CH_assert(new_puncture_coords.size() == m_num_punctures);

        for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
        {
            interp_x[ipuncture] = new_puncture_coords[ipuncture][0];
            interp_y[ipuncture] = new_puncture_coords[ipuncture][1];
            interp_z[ipuncture] = new_puncture_coords[ipuncture][2];
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

        // update the puncture tracks
        for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
        {
            new_puncture_coords[ipuncture][0] +=
                -m_dt * interp_shift1[ipuncture];
            new_puncture_coords[ipuncture][1] +=
                -m_dt * interp_shift2[ipuncture];
            new_puncture_coords[ipuncture][2] +=
                -m_dt * interp_shift3[ipuncture];
        }
        a_gramr.set_puncture_coords(new_puncture_coords);

        // print them out
        if (write_punctures)
        {
            SmallDataIO punctures_file(m_punctures_filename, m_dt, m_time,
                                       m_restart_time, SmallDataIO::APPEND);
            punctures_file.remove_duplicate_time_data();
            if (m_time == a_coarse_dt)
            {
                std::vector<std::string> header1_strings(CH_SPACEDIM *
                                                         m_num_punctures);
                header1_strings[0] = "x1";
                header1_strings[1] = "y1";
                header1_strings[2] = "z1";
                header1_strings[3] = "x2";
                header1_strings[4] = "y2";
                header1_strings[5] = "z2";
                punctures_file.write_header_line(header1_strings);
            }

            // use a vector for the write out
            std::vector<double> puncture_vector;
            puncture_vector.resize(m_num_punctures * CH_SPACEDIM);
            for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
            {
                puncture_vector[ipuncture * CH_SPACEDIM + 0] =
                    new_puncture_coords[ipuncture][0];
                puncture_vector[ipuncture * CH_SPACEDIM + 1] =
                    new_puncture_coords[ipuncture][1];
                puncture_vector[ipuncture * CH_SPACEDIM + 2] =
                    new_puncture_coords[ipuncture][2];
            }
            punctures_file.write_time_data_line(puncture_vector);
        }
    }
};

#endif /* PUNCTURETRACKER_HPP_ */
