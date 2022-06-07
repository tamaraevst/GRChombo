#ifndef COMPUTEWEIGHTFUNCTION_HPP_
#define COMPUTEWEIGHTFUNCTION_HPP_

// Chombo includes
#include "IntVect.H"

// #include "simd.hpp"
#include "VarsTools.hpp"
#include <array>
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DebuggingTools.hpp"
#include "BosonStarParams.hpp"
#include "WeightFunction.hpp"
#include "simd.hpp"

class ComputeWeightFunction
{
    protected:

    BosonStar_params_t m_params_BosonStar;
    BosonStar_params_t m_params_BosonStar2;
    double m_dx;

    public:

    template <class data_t> struct weightfunc_t
    {
        data_t weight1;
        data_t weight2;
    };

    ComputeWeightFunction(BosonStar_params_t a_params_BosonStar,
                   BosonStar_params_t a_params_BosonStar2, 
                   const double a_dx)
                    : m_dx(a_dx), 
                    m_params_BosonStar(a_params_BosonStar), 
                    m_params_BosonStar2(a_params_BosonStar2)
    {   
    }

    //! Function to compute the value of the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        weightfunc_t<data_t> out;

        Coordinates<double> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);

        double separation = m_params_BosonStar.BS_separation;
        double q = m_params_BosonStar.mass_ratio;
        double rapidity = m_params_BosonStar.BS_rapidity;
        double rapidity2 = m_params_BosonStar2.BS_rapidity;
        double impact_parameter = m_params_BosonStar.BS_impact_parameter;
        double radius_width1 = m_params_BosonStar.radius_width1;
        double radius_width2 = m_params_BosonStar.radius_width2;

        WeightFunction weightfunction;

        double profile_func1 = weightfunction.profile_chi((coords.x-separation/(q+1))*cosh(rapidity), coords.y, coords.z, radius_width1);
        double profile_func2 = weightfunction.profile_chi((coords.x+q*separation/(q+1))*cosh(-rapidity2), coords.y, coords.z, radius_width2);

        current_cell.store_vars(profile_func1, c_profile1);
        current_cell.store_vars(profile_func2, c_profile2);
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */
