/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHERNSIMONSEXTRACTION_HPP
#define CHERNSIMONSEXTRACTION_HPP

#include "SphericalExtraction.hpp"

//!  The class extracts the MatterEnergyFlux and integrates over spheres.
class ChernSimonsExtraction : public SphericalExtraction
{
  protected:
    int m_var_enum;
    VariableType m_var_type;

  public:
    //! The constructor
    ChernSimonsExtraction(
        SphericalExtraction::params_t &a_params, double a_dt, double a_time,
        bool a_first_step, double a_restart_time, int a_var_enum,
        VariableType a_var_type = VariableType::diagnostic)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time),
          m_var_enum(a_var_enum), m_var_type(a_var_type)
    {
        add_var(m_var_enum, m_var_type);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("ChernSimonsExtraction");
        }
    }
};

#endif /* CHERNSIMONSEXTRACTION_HPP */
