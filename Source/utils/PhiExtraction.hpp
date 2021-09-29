/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHIEXTRACTION_HPP_
#define PHIEXTRACTION_HPP_

#include "SphericalExtraction.hpp"

//!  The class allows extraction of the values of phi components on
//!  spherical shells at specified radii, and integration over those shells

class PhiExtraction : public SphericalExtraction
{
    // protected:
    // int m_var_enum;
    // VariableType m_var_type;

  public:
    //! The constructor
    // PhiExtraction(
    //     SphericalExtraction::params_t &a_params, double a_dt, double a_time,
    //     bool a_first_step, double a_restart_time = 0.0, int a_var_enum,
    //     VariableType a_var_type = VariableType::evolution)
    //     : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
    //                           a_restart_time),
    //       m_var_enum(a_var_enum), m_var_type(a_var_type)
    // {
    //     add_evolution_vars(m_var_enum, m_var_type);
    // }

    //! The constructor
    PhiExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_phi, VariableType::evolution);
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(m_params.extraction_file_prefix);

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<double>, std::vector<double>>>
            mode_integrals(m_num_modes);

         // note that this is normalised by multiplying by radius
        auto normalised_phi = [](std::vector<double> phi_reim_parts,
                                           double r, double, double) {
            // here the std::vector<double> passed will just have
            // the real and imaginary parts of the Weyl4 scalar as its
            // only components
            return std::make_pair(r * phi_reim_parts[0],
                                  0.0 * phi_reim_parts[1]);
        };

        // add the modes that will be integrated
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            constexpr int es = 0;
            add_mode_integrand(es, mode.first, mode.second,
                               normalised_phi, mode_integrals[imode]);
        }

        // do the integration over the surface
        integrate();

        // write the integrals
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            std::string integrals_phi_filename = m_params.integral_file_prefix +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_phi_for_writing = {
                std::move(mode_integrals[imode].first),
                std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral"};
            write_integrals(integrals_phi_filename, integrals_phi_for_writing, labels);
        }
    }
};

#endif /* PHIEXTRACTION_HPP_ */
