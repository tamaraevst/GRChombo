/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GAUSSIANFITTRACKING_HPP_
#define GAUSSIANFITTRACKING_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "UserVariables.hpp" // Needs c_mod_phi etc
#include "SimulationParametersBase.hpp"

//////
#include <fstream>
#include <string>
//////
class GaussianFitTracking
{
    //add linking to parameter file for search width and shit
    int m_number;
    int m_field_index;
    double m_delta, m_L;
    double* m_vals;
    double* m_array_x;
    double* m_array_y;
    double* m_array_z;
    double* m_centre = new double[3];
    double* m_old_centre = new double[3];
    double m_dt, m_time, m_restart_time;
    bool m_first_step;
    string m_filename = "StarTracking";
    std::ofstream m_test_file;
    GaussFit_params_t m_params_GaussFit;

  public:


    GaussianFitTracking(GaussFit_params_t a_params_GaussFit, double a_dt, double a_time,
                                    double a_restart_time, bool a_first_step, double a_L)
                                    :m_params_GaussFit(a_params_GaussFit),m_dt(a_dt),
                                    m_time(a_time),m_restart_time(a_restart_time),
                                    m_first_step(a_first_step),m_L(a_L)
    {
        m_number = m_params_GaussFit.num_points;
        m_field_index = m_params_GaussFit.field_index;
        m_delta = m_params_GaussFit.search_width;
        //m_test_file.open("TestDatFile.dat",std::ifstream::app);
        m_vals = new double[3*m_number];
        m_array_x = new double[3*m_number];
        m_array_y = new double[3*m_number];
        m_array_z = new double[3*m_number];
    }
    ~GaussianFitTracking()
    {
        delete[] m_vals;
        delete[] m_old_centre;
        delete[] m_centre;
        delete[] m_array_x;
        delete[] m_array_y;
        delete[] m_array_z;
        m_test_file.close();
    }
    void get_data(AMRInterpolator<Lagrange<4>> *a_interpolator);
    void find_centre();
    void write_to_dat();
    void do_star_tracking(AMRInterpolator<Lagrange<4>> *a_interpolator);
    void read_old_centre_from_dat();
};

#include "GaussianFitTracking.impl.hpp"

#endif /* GAUSSIANFITTRACKING_HPP_ */
