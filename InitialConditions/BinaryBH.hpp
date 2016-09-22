#ifndef BINARYBH_HPP_
#define BINARYBH_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "BoostedBH.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

namespace Lapse
{
    enum
    {
        ONE,
        PRE_COLLAPSED,
        CHI
    };
}

class BinaryBH
{
protected:
    const FABDriverBase& m_driver;
    double m_dx;
    BoostedBH bh1;
    BoostedBH bh2;
    int m_initial_lapse;

public:
    BinaryBH(const FABDriverBase& a_driver, BoostedBH::params_t a_bh1_params,
             BoostedBH::params_t a_bh2_params, double a_dx, int a_initial_lapse = Lapse::PRE_COLLAPSED) :
        m_driver (a_driver), m_dx (a_dx), bh1 (a_bh1_params), bh2 (a_bh2_params),
        m_initial_lapse (a_initial_lapse) {}

    //Not currently vectorised (it is only done once so it's hardly worth adding all the special functions to simd)
    void compute(int ix, int iy, int iz);

protected:
    double compute_chi(double x, double y, double z);

    tensor<2,double> compute_A(double chi, double x, double y, double z);

};

void BinaryBH::compute(int ix, int iy, int iz)
{
    CCZ4::vars_t<double> vars;
    vars.assign(0.); //Set only the non-zero components explicitly below
    Coordinates<double> coords(ix,iy,iz,m_dx);
    double x = coords.x; //TODO: change functions to accept coords rather than x,y,z
    double y = coords.y;
    double z = coords.z;

    vars.chi = compute_chi(x,y,z);

    //Conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    vars.A = compute_A(vars.chi,x,y,z);

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

    idx_t<double> out_idx = m_driver.out_idx(ix, iy, iz); //The current location in the flattened output FArraBox
    m_driver.store_vars(vars, out_idx);
}

double BinaryBH::compute_chi(double x, double y, double z)
{
    const double psi = 1. + bh1.psi_minus_one(x, y, z) + bh2.psi_minus_one(x, y, z);
    return pow(psi, -4);
}

tensor<2,double> BinaryBH::compute_A(double chi, double x, double y, double z)
{

    tensor<2,double> Aij1 = bh1.Aij(x,y,z);
    tensor<2,double> Aij2 = bh2.Aij(x,y,z);
    tensor<2, double> out;

    //Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR2(i,j) out[i][j] = pow(chi, 3/2.) * (Aij1[i][j] + Aij2[i][j]);

    return out;
}

#endif /* BINARYBH_HPP_ */
