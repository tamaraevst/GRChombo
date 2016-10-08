#ifndef BUBBLESF_HPP_
#define BUBBLESF_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "FABDriverBase.hpp"
#include <vector>
#include "tensor.hpp"
#include <array>
#include "UserVariables.hpp" //This files needs c_NUM - total number of components

class BubbleSF
{
protected:
    const FABDriverBase& m_driver;
    double m_dx;

public:
    struct params_t
    {
        double amplitudeSF;
        std::vector<double> centerSF;
        double widthSF;
    };

    const params_t m_params;

    BubbleSF(const FABDriverBase& a_driver, params_t a_params,
             double a_dx) :
        m_driver (a_driver), m_dx (a_dx), m_params (a_params) {}

    //Not currently vectorised (it is only done once so it's hardly worth adding all the special functions to simd)
    void compute(int ix, int iy, int iz);

protected:
    double compute_phi(double x, double y, double z);

};

#endif /* BUBBLESF_HPP_ */
