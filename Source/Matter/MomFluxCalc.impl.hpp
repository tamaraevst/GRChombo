/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MOMFLUXCALC_HPP)
#error "This file should only be included through MomFluxCalc.hpp"
#endif

#ifndef MOMFLUXCALC_IMPL_HPP
#define MOMFLUXCALC_IMPL_HPP

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "simd.hpp"


template <class matter_t>
EMTensor_and_mom_flux<matter_t>::EMTensor_and_mom_flux(const matter_t &a_matter,
                             const double dx, const double a_L,
                              std::array<double,CH_SPACEDIM> a_centre, const int a_c_Fx_flux,
                             const int a_c_Fy_flux, const int a_c_rho,
                             const Interval a_c_Si, const Interval a_c_Sij)
    : m_matter(a_matter), m_deriv(dx),  m_dx(dx), m_L(a_L),
                                  m_centre(a_centre), m_c_Fx_flux(a_c_Fx_flux),
                                   m_c_Fy_flux(a_c_Fy_flux),  m_c_rho(a_c_rho),
                                              m_c_Si(a_c_Si), m_c_Sij(a_c_Sij)
{
    if (m_c_Si.size() != 0)
    {
        // Si is a vector
        CH_assert(m_c_Si.size() == DEFAULT_TENSOR_DIM);
    }

    if (m_c_Sij.size() != 0)
    {
        // Sij is a symmetric tensor
        CH_assert(m_c_Sij.size() ==
                  DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2);
    }
}

template <class matter_t>
template <class data_t>
void EMTensor_and_mom_flux<matter_t>::compute(Cell<data_t> current_cell) const
{
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    if (m_c_rho >= 0)
    {
        current_cell.store_vars(emtensor.rho, m_c_rho);
    }

    if (m_c_Si.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        FOR1(i) { current_cell.store_vars(emtensor.Si[i], m_c_Si.begin() + i); }
#endif
    }

    if (m_c_Sij.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        current_cell.store_vars(emtensor.Sij[0][0], m_c_Sij.begin());
        current_cell.store_vars(emtensor.Sij[0][1], m_c_Sij.begin() + 1);
        current_cell.store_vars(emtensor.Sij[0][2], m_c_Sij.begin() + 2);
        current_cell.store_vars(emtensor.Sij[1][1], m_c_Sij.begin() + 3);
        current_cell.store_vars(emtensor.Sij[1][2], m_c_Sij.begin() + 4);
        current_cell.store_vars(emtensor.Sij[2][2], m_c_Sij.begin() + 5);
#endif
    }

    Coordinates<data_t> coords(current_cell, m_dx,m_centre);

    //double N[3] = {coords.x,coords.y,coords.z};
    data_t N[3] = {coords.x,coords.y,coords.z};
    data_t normsqr = 0.;
    data_t Fx=0., Fy=0.;

    FOR2(i,j) normsqr += N[i]*N[j]*vars.h[i][j]/vars.chi;
    FOR1(i) N[i] /= sqrt(normsqr);
    FOR1(i) Fx += vars.lapse*N[i]*emtensor.Sij[i][0];
    FOR1(i) Fy += vars.lapse*N[i]*emtensor.Sij[i][1];
    FOR2(i,j) Fx -=  vars.h[i][j]*vars.shift[i]*N[j]*emtensor.Si[0]/vars.chi;
    FOR2(i,j) Fy -=  vars.h[i][j]*vars.shift[i]*N[j]*emtensor.Si[1]/vars.chi;


    current_cell.store_vars(Fx,m_c_Fx_flux);
    current_cell.store_vars(Fy,m_c_Fy_flux);

}

#endif /* MOMFLUXCALC_IMPL_HPP */
