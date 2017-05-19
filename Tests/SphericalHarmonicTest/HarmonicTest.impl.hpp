#if !defined(HARMONICTEST_HPP_)
#error "This file should only be included through HarmonicTest.hpp"
#endif

#ifndef HARMONICTEST_IMPL_HPP_
#define HARMONICTEST_IMPL_HPP_

#include "HarmonicTest.hpp"
#include "ScalarField.hpp"
#include "SphericalHarmonics.hpp"
#include "simd.hpp"
#include "DebuggingTools.hpp"

template <class data_t>
void HarmonicTest::compute(Cell current_cell) {

    ScalarField<>::Vars<data_t> vars;
    Coordinates<data_t> coords(current_cell, m_dx);

    vars.phi = compute_harmonic(coords);

    data_t radius = Coordinates<data_t>::get_radius(current_cell, m_dx, m_center_vector);
    vars.phi = vars.phi/radius;

    current_cell.store_vars(vars.phi, c_phi);
}

template <class data_t>
data_t HarmonicTest::compute_harmonic(Coordinates<data_t> coords) {

    // work out where we are on the grid
    data_t x = coords.x - m_center_vector[0];
    double y = coords.y - m_center_vector[1];
    double z = coords.z - m_center_vector[2];

    data_t rr = coords.get_radius(m_center_vector);

    //Add in el, em spherical harmonics here, spin weight es
    using namespace SphericalHarmonics;
    int es = -1;
    int el = 2;
    int em = -1;
    auto Y_lm = spin_Y_lm(x, y, z, rr, es, el, em);
    data_t out = Y_lm.Real;

    return out;
}

#endif /* HARMONICTEST_IMPL_HPP_ */

