#ifndef BINARYBH_HPP_
#define BINARYBH_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "BoostedBH.hpp"
#include "Cell.hpp"

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
    double m_dx;
    BoostedBH bh1;
    BoostedBH bh2;
    int m_initial_lapse;

public:
    BinaryBH(BoostedBH::params_t a_bh1_params, BoostedBH::params_t a_bh2_params, double a_dx,
             int a_initial_lapse = Lapse::PRE_COLLAPSED) :
        m_dx (a_dx), bh1 (a_bh1_params), bh2 (a_bh2_params),
        m_initial_lapse (a_initial_lapse) {}

    template <class data_t>
    void compute(Cell current_cell);

protected:
    template <class data_t>
    data_t compute_chi(Coordinates<data_t> coords);

    template <class data_t>
    tensor<2,data_t> compute_A(data_t chi, Coordinates<data_t> coords);

};

#include "BinaryBH.impl.hpp"

#endif /* BINARYBH_HPP_ */
