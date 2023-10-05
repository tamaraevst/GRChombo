/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_
#define BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

//! This class tags cells based on three criteria - the
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to
//! a given level
class BosonChiPunctureExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_max_level;
    const bool m_track_punctures;
    const bool m_activate_extraction;
    const FourthOrderDerivatives m_deriv;
    const SphericalExtraction::params_t m_params;
    const std::vector<double> m_puncture_radii;
    const std::vector<double> m_puncture_masses;
    const std::vector<std::array<double, CH_SPACEDIM>> &m_puncture_coords;
    const double m_buffer;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    BosonChiPunctureExtractionTaggingCriterion(
        const double dx, const int a_level, const int a_max_level,
        const SphericalExtraction::params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> &a_puncture_coords,
        const bool activate_extraction = false,
        const bool track_punctures = false,
        const std::vector<double> a_puncture_radii = {4.0, 4.0}, 
        const std::vector<double> a_puncture_masses = {1.0, 1.0}, 
        const double a_buffer = 0.5)
        : m_dx(dx), m_level(a_level), m_max_level(a_max_level),
          m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction), m_deriv(dx),
          m_params(a_params), m_puncture_radii(a_puncture_radii), 
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords),
          m_buffer(a_buffer)
    {
        // check that the number of punctures is consistent
        CH_assert(m_puncture_radii.size() == m_puncture_coords.size());
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        if (m_track_punctures == 1)
        {
            // we want the 2nd and 3rd finest levels to be twice the size of the
            // next finest level
            const double factor = pow(2.0, min(m_max_level - m_level - 1, 2));
            // loop over puncture numbers
            for (int ipuncture = 0; ipuncture < m_puncture_radii.size(); ++ipuncture)
            {
                // where am i?
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_puncture_coords[ipuncture]);
                const data_t r = coords.get_radius();

                auto regrid =
                    simd_compare_lt(r, (1.0 + m_buffer) * factor *
                                           m_puncture_masses[ipuncture]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }

            double puncture_separation2 = 0.0;
            FOR1(idir)
                {
                    double displacement = m_puncture_coords[0][idir] - m_puncture_coords[1][idir];
                    puncture_separation2 += displacement * displacement;
                }

            double puncture_separation = sqrt(puncture_separation2);

            double sum_radii = m_puncture_radii[0] + m_puncture_radii[1];
            double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];

            // if punctures are close enough together tag cells at the
            // center of mass for the merger BH
            if (puncture_separation < (1.0 + m_buffer) * sum_masses)
            {
                std::array<double, CH_SPACEDIM> center_of_mass;
                FOR1(idir)
                {
                    center_of_mass[idir] =
                    (m_puncture_masses[0] * m_puncture_coords[0][idir] +
                    m_puncture_masses[1] * m_puncture_coords[1][idir]) / sum_masses;
                }
            
                Coordinates<data_t> coords(current_cell, m_dx, center_of_mass);
                const data_t r = coords.get_radius();

                auto regrid2 = simd_compare_lt(r, (1.0 + m_buffer) * factor * sum_masses);
                criterion = simd_conditional(regrid2, 100.0, criterion);
            }
        }

        // ensure that the horizons of the punctures are covered
        // by the max level - for this we need
        // only check the puncture locations on the top 2 levels
        // which regrid (ie, max_level - 1 to max_level - 2)
        // (just the top level would be ok, but doing two ensures
        // the top levels are well spaced)
        if ((m_level > (m_max_level - 3)) && (m_track_punctures == 1))
        {
            // we want each level to be double the innermost one in size
            const double factor = pow(2.0, m_max_level - m_level - 1);
            // loop over puncture radii
            for (int ipuncture = 0; ipuncture < m_puncture_radii.size();
                 ++ipuncture)
            {
                // where am i?
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_puncture_coords[ipuncture]);
                const data_t r = coords.get_radius();
                // decide whether to tag based on distance to horizon
                // plus a fudge factor of 1.5
                auto regrid = simd_compare_lt(
                    r, 1.5 * factor * m_puncture_radii[ipuncture]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BOSONCHIPUNCTUREEXTRACTIONTAGGINGCRITERION_HPP_ */
