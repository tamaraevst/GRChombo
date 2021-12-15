#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

// Chombo includes
#include "IntVect.H"

#include "simd.hpp"
#include "VarsTools.hpp"
#include <array>
#include "Cell.hpp"
#include "Coordinates.hpp"

class WeightFunction
{
    protected:
    double m_dx;
    std::array<double, CH_SPACEDIM> m_center;

    public:
    WeightFunction(const double a_dx, const std::array<double, CH_SPACEDIM> a_center) : m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t>
    double compute_weight(Cell<data_t> current_cell, double pos, double width) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        
        double weightfunc;

        if (coords.z = 0.0)
        {
            double scaledr = sqrt(pow(coords.x - pos,2) + pow(coords.y, 2)) / (width) ;

            if (scaledr <= 1)
            {
                return weightfunc = 140 * ( (1.0/4.0)*pow((1-scaledr),4) - (3.0/5.0)*pow((1-scaledr),5) + (1.0/2.0)* pow((1-scaledr),6) - (1.0)/(7.0)*pow((1-scaledr),7));
            }

            if (scaledr > 1)
            {
                return weightfunc = 0.0;
            }
        }
        else
        {
            return weightfunc = 0.0;
        }
    }
    
};

#endif /* WEIGHTFUNCTION_HPP_ */