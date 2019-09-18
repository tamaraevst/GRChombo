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
    //! Params for extraction
    const int m_num_punctures;
    const double m_time;
    const double m_restart_time;
    const double m_dt;
    const std::string m_checkpoint_prefix;

  public:
    //! The constructor
    PunctureTracker(const double a_time, const double a_restart_time,
                    const double a_dt, const std::string a_checkpoint_prefix,
                    const int a_num_punctures = 2)
        : m_num_punctures(a_num_punctures), m_time(a_time),
          m_checkpoint_prefix(a_checkpoint_prefix),
          m_restart_time(a_restart_time), m_dt(a_dt)
    {
    }

    //! Execute the query
    void execute_tracking(GRAMR &a_gramr) const
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
        std::vector<double> new_puncture_coords = a_gramr.get_puncture_coords();
        CH_assert(new_puncture_coords.size() / CH_SPACEDIM == m_num_punctures);
        for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
        {
            interp_x[ipuncture] = new_puncture_coords[CH_SPACEDIM * ipuncture];
            interp_y[ipuncture] =
                new_puncture_coords[CH_SPACEDIM * ipuncture + 1];
            interp_z[ipuncture] =
                new_puncture_coords[CH_SPACEDIM * ipuncture + 2];
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
            new_puncture_coords[CH_SPACEDIM * ipuncture] +=
                -m_dt * interp_shift1[ipuncture];
            new_puncture_coords[CH_SPACEDIM * ipuncture + 1] +=
                -m_dt * interp_shift2[ipuncture];
            new_puncture_coords[CH_SPACEDIM * ipuncture + 2] +=
                -m_dt * interp_shift3[ipuncture];
        }
        a_gramr.set_puncture_coords(new_puncture_coords);

        // print them out
        std::string extraction_filename =
            m_checkpoint_prefix + "_Punctures.txt";
        SmallDataIO extraction_file(extraction_filename, m_dt, m_time,
                                    m_restart_time, SmallDataIO::APPEND);
        extraction_file.remove_duplicate_time_data();
        if (m_time == m_dt)
        {
            std::vector<std::string> header1_strings(CH_SPACEDIM *
                                                     m_num_punctures);
            header1_strings[0] = "x1";
            header1_strings[1] = "y1";
            header1_strings[2] = "z1";
            header1_strings[3] = "x2";
            header1_strings[4] = "y2";
            header1_strings[5] = "z2";
            extraction_file.write_header_line(header1_strings);
        }
        extraction_file.write_time_data_line(new_puncture_coords);
    }
};

#endif /* PUNCTURETRACKER_HPP_ */
