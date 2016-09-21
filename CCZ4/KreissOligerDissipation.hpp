#ifndef KREISSOLIGERDISSIPATION_HPP_
#define KREISSOLIGERDISSIPATION_HPP_

#include "simd.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

class KreissOligerDissipation
{
protected:
    const double m_sigma;
    const FABDriverBase& m_driver;
    const FourthOrderDerivatives m_deriv;

public:
    KreissOligerDissipation(double sigma, double dx, const FABDriverBase& driver) :
        m_sigma (sigma),
        m_driver (driver),
        m_deriv (dx, m_driver)
    {}

    template <class data_t>
    void compute(int ix, int iy, int iz)
    {
        data_t sigma = m_sigma;

        idx_t<data_t> in_idx = m_driver.in_idx(ix, iy, iz);
        auto dissipation = m_deriv.dissipation(in_idx);

        idx_t<data_t> out_idx = m_driver.out_idx(ix, iy, iz);
        for (int comp = 0; comp < c_NUM; ++comp)
        {
            SIMDIFY<data_t>(m_driver.m_out_ptr[comp])[out_idx] += sigma * dissipation[comp];
        }
    }
};

#endif /* KREISSOLIGERDISSIPATION_HPP_ */
