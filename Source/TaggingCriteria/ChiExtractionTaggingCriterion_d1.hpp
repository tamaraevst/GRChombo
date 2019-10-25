/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIEXTRACTIONTAGGINGD1CRITERION_HPP_
#define CHIEXTRACTIONTAGGINGD1CRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

//! This class tags cells based on three criteria - the 
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to 
//! a given level
class ChiExtractionTaggingCriterion_d1
{
  protected:
    const double m_dx;
    const int m_level;
    const bool m_activate_extraction;
    const extraction_params_t m_params;
    const FourthOrderDerivatives m_deriv;

  public:
    // The constructor
    ChiExtractionTaggingCriterion_d1(const double dx, const int a_level,
                                  const extraction_params_t a_params,
                                  const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx), m_params(a_params), m_level(a_level),
          m_activate_extraction(activate_extraction)
    {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto chi = current_cell.load_vars(c_chi);
        Tensor<1, data_t> d1_chi;
        FOR1(idir) m_deriv.diff1(d1_chi, current_cell, idir, c_chi);

        data_t mod_d1_chi = 0;
        FOR1(idir) mod_d1_chi += d1_chi[idir] * d1_chi[idir];
        data_t criterion = m_dx * sqrt(mod_d1_chi) / pow(chi, 2);

        // if extracting weyl data at a given radius, enforce a given resolution there
        if (m_activate_extraction)
        { 
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
        }
    
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIEXTRACTIONTAGGINGD1CRITERION_HPP_ */
