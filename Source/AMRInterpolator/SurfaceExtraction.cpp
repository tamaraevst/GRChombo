/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SurfaceExtraction.hpp"

//! Normal constructor which requires vars to be added after construction
//! using add_var or add_vars
SurfaceExtraction::SurfaceExtraction(SurfaceGeometry *a_geom_ptr,
                                     params_t &a_params, double a_dt,
                                     double a_time, bool a_first_step,
                                     double a_restart_time)
    : m_geom_ptr(a_geom_ptr), m_params(a_params), m_dt(a_dt), m_time(a_time),
      m_first_step(a_first_step), m_restart_time(a_restart_time),
      m_num_points(m_params.num_points_u * m_params.num_points_v),
      m_du(m_geom_ptr->du(m_params.num_points_u)),
      m_dv(m_geom_ptr->dv(m_params.num_points_v)), m_done_extraction(false)
{
    if (m_geom_ptr == nullptr)
    {
        MayDay::Error("SurfaceExtraction: invalid geometry pointer");
    }
    FOR1(idir)
    {
        m_interp_coords[idir].resize(m_num_points * m_params.num_surfaces);
    }

    for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
    {
        double surface_param_value = m_params.surface_param_values[isurface];
        for (int iu = 0; iu < m_params.num_points_u; ++iu)
        {
            double u = m_geom_ptr->u(iu, m_params.num_points_u);
            for (int iv = 0; iv < m_params.num_points_v; ++iv)
            {
                double v = m_geom_ptr->v(iv, m_params.num_points_v);
                FOR1(idir)
                {
                    int idx = index(isurface, iu, iv);
                    m_interp_coords[idir][idx] = m_geom_ptr->cartesian_coord(
                        idir, surface_param_value, u, v);
                }
            }
        }
    }
}

//! add a single variable or derivative of variable
void SurfaceExtraction::add_var(const int a_var, const Derivative a_deriv)
{
    CH_assert(!m_done_extraction);
    m_vars.push_back({a_var, a_deriv});
    m_interp_data.emplace_back(m_num_points * m_params.num_surfaces);
}

//! add a vector of variables/derivatives of variables
void SurfaceExtraction::add_vars(
    const std::vector<std::pair<int, Derivative>> &a_vars)
{
    for (auto var : a_vars)
    {
        add_var(var.first, var.second);
    }
}

//! Alternative constructor with a predefined vector of variables
SurfaceExtraction::SurfaceExtraction(
    SurfaceGeometry *a_geom_ptr, params_t &a_params,
    const std::vector<std::pair<int, Derivative>> &a_vars, double a_dt,
    double a_time, bool a_first_step, double a_restart_time)
    : SurfaceExtraction(a_geom_ptr, a_params, a_dt, a_time, a_first_step,
                        a_restart_time)
{
    add_vars(a_vars);
}

// define the static IntegrationMethods here
const IntegrationMethod IntegrationMethod::trapezium({0.5});
const IntegrationMethod IntegrationMethod::midpoint({1.0}, false);
const IntegrationMethod IntegrationMethod::simpson({0.3333333333333333,
                                                     1.3333333333333333});
const IntegrationMethod
    IntegrationMethod::boole({0.3111111111111111, 1.4222222222222222,
                              0.53333333333333, 1.4222222222222222});

//! Integrate some integrand dependent on the interpolated data over the
//! surface. The integrand function should be of the signature
//! double integrand(std::vector<double> data_here,
//!     double a_surface_param_value, double a_u, double a_v)
//! where data_here is a vector of all the interpolated variables at the
//! point specified by the other arguments.
std::vector<double> SurfaceExtraction::integrate(
    std::function<double(std::vector<double>, double, double, double)>
        a_integrand,
    const IntegrationMethod &a_method_u, const IntegrationMethod &a_method_v)
{
    CH_assert(m_done_extraction);
    std::vector<double> out_integrals(m_params.num_surfaces, 0.0);

    bool valid_u =
        a_method_u.is_valid(m_params.num_points_u, m_geom_ptr->is_u_periodic());
    bool valid_v =
        a_method_v.is_valid(m_params.num_points_v, m_geom_ptr->is_v_periodic());

    // default to using the trapezium rule if provided methods are not valid
    IntegrationMethod method_u = IntegrationMethod::trapezium;
    IntegrationMethod method_v = IntegrationMethod::trapezium;
    if (!valid_u)
    {
        MayDay::Warning(
            "SurfaceExtraction::integrate: Provided "
            "IntegrationMethod for u is not valid with\nthis num_points_u; "
            "reverting to trapezium rule.");
    }
    else
    {
        method_u = a_method_u;
    }
    if (!valid_v)
    {
        MayDay::Warning(
            "SurfaceExtraction::integrate: Provided "
            "IntegrationMethod for v is not valid with\nthis num_points_v; "
            "reverting to trapezium rule.");
    }
    else
    {
        method_v = a_method_v;
    }

    for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
    {
        double surface_param_value = m_params.surface_param_values[isurface];
        for (int iu = 0; iu < m_params.num_points_u; ++iu)
        {
            double u = m_geom_ptr->u(iu, m_params.num_points_u);
            double inner_integral = 0.0;
            for (int iv = 0; iv < m_params.num_points_v; ++iv)
            {
                double v = m_geom_ptr->v(iv, m_params.num_points_v);
                std::vector<double> data_here(m_vars.size());
                for (int ivar = 0; ivar < m_vars.size(); ++ivar)
                {
                    data_here[ivar] =
                        m_interp_data[ivar][index(isurface, iu, iv)];
                }
                double integrand_with_area_element =
                    a_integrand(data_here, surface_param_value, u, v) *
                    m_geom_ptr->area_element(surface_param_value, u, v);
                double weight = method_v.weight(iv, m_params.num_points_v,
                                                m_geom_ptr->is_v_periodic());
                inner_integral += weight * m_dv * integrand_with_area_element;
            }
            double weight = method_u.weight(iu, m_params.num_points_u,
                                            m_geom_ptr->is_u_periodic());
            out_integrals[isurface] += weight * m_du * inner_integral;
        }
    }
    return out_integrals;
}

//! Write the interpolated data to a file with a block for each surface
void SurfaceExtraction::write_extraction(std::string a_file_prefix) const
{
    CH_assert(m_done_extraction);
    SmallDataIO extraction_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW, m_first_step);

    for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
    {
        // Write headers
        std::vector<std::string> header1_strings = {
            "time = " + std::to_string(m_time) + ",",
            m_geom_ptr->param_name() + " = " +
                std::to_string(m_params.surface_param_values[isurface])};
        extraction_file.write_header_line(header1_strings, "");
        std::vector<std::string> components(m_vars.size());
        for (int ivar = 0; ivar < m_vars.size(); ++ivar)
        {
            if (m_vars[ivar].second != Derivative::LOCAL)
            {
                components[ivar] = Derivative::name(m_vars[ivar].second) + "_";
            }
            else
            {
                components[ivar] = "";
            }
            components[ivar] +=
                UserVariables::variable_names[m_vars[ivar].first];
        }
        std::vector<std::string> coords = {m_geom_ptr->u_name(),
                                           m_geom_ptr->v_name()};
        extraction_file.write_header_line(components, coords);

        // Now the data
        for (int iu = 0; iu < m_params.num_points_u; ++iu)
        {
            double u = m_geom_ptr->u(iu, m_params.num_points_u);
            for (int iv = 0; iv < m_params.num_points_v; ++iv)
            {
                double v = m_geom_ptr->v(iv, m_params.num_points_v);
                int idx = index(isurface, iu, iv);
                std::vector<double> data(m_vars.size());
                for (int ivar = 0; ivar < m_vars.size(); ++ivar)
                {
                    data[ivar] = m_interp_data[ivar][idx];
                }

                extraction_file.write_data_line(data, {u, v});
            }
        }
        extraction_file.line_break();
    }
}

//! write some integrals to a file at this timestep
void SurfaceExtraction::write_integrals(
    const std::string &a_filename,
    const std::vector<std::vector<double>> &a_integrals,
    const std::vector<std::string> &a_labels) const
{
    const int num_integrals_per_surface = a_integrals.size();
    // if labels are provided there must be the same number of labels as
    // there are integrals
    if (!a_labels.empty())
    {
        CH_assert(num_integrals_per_surface == a_labels.size());
    }
    // each inner vector element of a_integrals must have the same number of
    // elements as there are surfaces (i.e. one integral per surface)
    for (auto vect : a_integrals)
    {
        CH_assert(vect.size() == m_params.num_surfaces);
    }
    // open file for writing
    SmallDataIO integral_file(a_filename, m_dt, m_time, m_restart_time,
                              SmallDataIO::APPEND, m_first_step);

    // remove any duplicate data if this is a restart
    integral_file.remove_duplicate_time_data();

    if (m_first_step)
    {
        // make header strings
        std::vector<std::string> header1_strings(num_integrals_per_surface *
                                                 m_params.num_surfaces);
        std::vector<std::string> header2_strings(num_integrals_per_surface *
                                                 m_params.num_surfaces);
        for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
        {
            for (int iintegral = 0; iintegral < num_integrals_per_surface;
                 ++iintegral)
            {
                int idx = isurface * num_integrals_per_surface + iintegral;
                if (a_labels.empty())
                    header1_strings[idx] = "";
                else
                    header1_strings[idx] = a_labels[iintegral];
                header2_strings[idx] =
                    std::to_string(m_params.surface_param_values[isurface]);
            }
        }
        std::string pre_header2_string = m_geom_ptr->param_name() + " = ";

        // write headers
        integral_file.write_header_line(header1_strings);
        integral_file.write_header_line(header2_strings, pre_header2_string);
    }

    // make vector of data for writing
    std::vector<double> data_for_writing(num_integrals_per_surface *
                                         m_params.num_surfaces);
    for (int isurface = 0; isurface < m_params.num_surfaces; ++isurface)
    {
        for (int iintegral = 0; iintegral < num_integrals_per_surface;
             ++iintegral)
        {
            int idx = isurface * num_integrals_per_surface + iintegral;
            data_for_writing[idx] = a_integrals[iintegral][isurface];
        }
    }

    // write data
    integral_file.write_time_data_line(data_for_writing);
}

//! convenience caller for write_integrals in the case of just one integral per
//! surface
void SurfaceExtraction::write_integral(const std::string &a_filename,
                                       const std::vector<double> a_integrals,
                                       const std::string a_label) const
{
    std::vector<std::vector<double>> integrals(1, a_integrals);
    if (!a_label.empty())
    {
        std::vector<std::string> labels(1, a_label);
        write_integrals(a_filename, integrals, labels);
    }
    else
        write_integrals(a_filename, integrals);
}
