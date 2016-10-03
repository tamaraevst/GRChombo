#ifndef BINARYBH_HPP_
#define BINARYBH_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "BoostedBH.hpp"
#include "FABDriverBase.hpp"

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

#endif /* BINARYBH_HPP_ */
