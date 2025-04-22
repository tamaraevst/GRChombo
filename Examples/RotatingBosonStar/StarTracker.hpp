/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STARTRACKER_HPP_
#define STARTRACKER_HPP_

#include "AMRInterpolator.hpp"
#include "AlwaysInline.hpp"
#include "Lagrange.hpp"

//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class StarTracker
{
  private:
    //! Params for puncture tracking
    std::array<double, CH_SPACEDIM> m_star_coords;
    std::array<double, CH_SPACEDIM> m_centre;
    int m_tracking_level; // level (i.e. times) to execute tracking
    int m_points;     // number of points n, (2n + 1 points to integrate)
    double m_width_A;
    std::string m_direction;

    // saved pointer to external interpolator
    AMRInterpolator<Lagrange<4>> *m_interpolator;

  public:
    //! The constructor
    StarTracker() : m_interpolator(nullptr) {}

    //! set puncture locations on start (or restart)
    //! this needs to be done before 'setupAMRObject'
    //! if the puncture locations are required for Tagging Criteria
    void initial_setup(bool a_do_star_track, std::array<double, CH_SPACEDIM> &a_initial_star_centres,
                       int a_star_points, double a_star_track_width_A, std::string a_direction)
    {	
        
      m_points = a_star_points;
      m_width_A = a_star_track_width_A;
      m_direction = a_direction;

      for (int i = 0; i < CH_SPACEDIM; i++)
          {
            m_star_coords[i] = a_initial_star_centres[i];
            pout() << "\n Initialising the coordinate number " << i << " at " << a_initial_star_centres[i] << "\n" << std::endl;
          }
    }

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    double find_centre(int direction);

    void find_max_min(int direction);

    void update_star_centres(double a_dt);

    void write_to_dat(std::string a_filename, double a_dt, double a_time,
                      double a_restart_time, bool a_first_step);

    void read_old_centre_from_dat(std::string a_filename, double a_dt,
                                           double a_time, double a_restart_time,
                                           bool a_first_step);
    // function to get punctures
    ALWAYS_INLINE const std::array<double, CH_SPACEDIM> &
    get_puncture_coords() const
    {
        return m_star_coords;
    }
};

#endif /* STARTRACKER_HPP_ */
