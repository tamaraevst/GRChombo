#include "BinaryBH.hpp"
#include "CCZ4.hpp"

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

    m_driver.store_vars(vars);
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
