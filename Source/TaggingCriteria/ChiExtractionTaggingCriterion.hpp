/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIEXTRACTIONTAGGINGCRITERION_HPP_
#define CHIEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class ChiExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_max_level;
    const int m_track_punctures;
    const int m_activate_extraction;
    const extraction_params_t m_params;
    const FourthOrderDerivatives m_deriv;
    std::array<double, 3> m_puncture_coords1;
    std::array<double, 3> m_puncture_coords2;

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

    ChiExtractionTaggingCriterion(double dx, int a_level, int a_max_level,
                                  extraction_params_t a_params,
                                  std::vector<double> a_puncture_coords,
                                  int activate_extraction = 0,
                                  int track_punctures = 0)
        : m_dx(dx), m_deriv(dx), m_params(a_params), m_level(a_level),
          m_max_level(a_max_level), m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction){
              FOR1(i){m_puncture_coords1[i] = a_puncture_coords[i];
    m_puncture_coords2[i] = a_puncture_coords[3 + i];
}
}
;

template <class data_t> void compute(Cell<data_t> current_cell) const
{
    const auto d2 = m_deriv.template diff2<Vars>(current_cell);

    data_t mod_d2_chi = 0;
    FOR2(idir, jdir) { mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir]; }

    data_t criterion = m_dx * sqrt(mod_d2_chi);

    for (int iradius = 0; iradius < m_params.num_extraction_radii; ++iradius)
    {
        // regrid if within extraction level and not at required refinement
        if (m_level < m_params.extraction_levels[iradius])
        {
            const Coordinates<data_t> coords(current_cell, m_dx,
                                             m_params.extraction_center);
            const data_t r = coords.get_radius();
            // add a 20% buffer to extraction zone so not too near to
            // boundary
            auto regrid =
                simd_compare_lt(r, 1.2 * m_params.extraction_radii[iradius]);
            criterion = simd_conditional(regrid, 100.0, criterion);
        }
    }

    // only check the puncture locations on the top level
    if ((m_level > (m_max_level - 1)) && (m_track_punctures == 1))
    {
        const Coordinates<data_t> coords1(current_cell, m_dx,
                                          m_puncture_coords1);
        const Coordinates<data_t> coords2(current_cell, m_dx,
                                          m_puncture_coords2);
        const data_t r1 = coords1.get_radius();
        const data_t r2 = coords2.get_radius();
        auto regrid = simd_compare_lt(r1, 0.75);
        criterion = simd_conditional(regrid, 1.0e10, criterion);
        regrid = simd_compare_lt(r2, 0.75);
        criterion = simd_conditional(regrid, 1.0e10, criterion);
    }

    // Write back into the flattened Chombo box
    current_cell.store_vars(criterion, 0);
}
}
;

#endif /* CHIEXTRACTIONTAGGINGCRITERION_HPP_ */
