#if !defined(BINARYBH_HPP_)
#error "This file should only be included through BinaryBH.hpp"
#endif

#ifndef BINARYBH_IMPL_HPP_
#define BINARYBH_IMPL_HPP_

#include "BinaryBH.hpp"
#include "CCZ4.hpp"
#include "simd.hpp"
#include "Cell.hpp"

template <class data_t>
void BinaryBH::compute(Cell current_cell)
{
    CCZ4::Vars<data_t> vars;
    vars.assign(0.); //Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell,m_dx);

    vars.chi = compute_chi(coords);

    //Conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    vars.A = compute_A(vars.chi,coords);

    switch (m_initial_lapse)
    {
        case Lapse::ONE:
            vars.lapse = 1.;
            break;
        case Lapse::PRE_COLLAPSED:
            vars.lapse = sqrt(vars.chi);
            break;
        case Lapse::CHI:
            vars.lapse = vars.chi;
            break;
        default:
            MayDay::Error("BinaryBH::Supplied initial lapse not supported.");
    }

    m_driver.store_vars(vars, current_cell);
}

template <class data_t>
data_t BinaryBH::compute_chi(Coordinates<data_t> coords)
{
    const data_t psi = 1. + bh1.psi_minus_one(coords) + bh2.psi_minus_one(coords);
    return pow(psi, -4);
}

template <class data_t>
tensor<2,data_t> BinaryBH::compute_A(data_t chi, Coordinates<data_t> coords)
{

    tensor<2,data_t> Aij1 = bh1.Aij(coords);
    tensor<2,data_t> Aij2 = bh2.Aij(coords);
    tensor<2, data_t> out;

    //Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR2(i,j) out[i][j] = pow(chi, 3/2.) * (Aij1[i][j] + Aij2[i][j]);

    return out;
}

#endif /* BINARYBH_IMPL_HPP_ */
