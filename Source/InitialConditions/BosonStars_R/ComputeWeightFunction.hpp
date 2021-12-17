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
    const std::array<double, CH_SPACEDIM> m_center;

    public:

    template <class data_t> struct weightfunc_t
    {
        data_t weight1;
        data_t weight2;
    };

    ComputeWeightFunction(BosonStar_params_t a_params_BosonStar,
                   BosonStar_params_t a_params_BosonStar2, 
                   const double a_dx, 
                   const std::array<double, CH_SPACEDIM> a_center)
                    : m_dx(a_dx), m_center(a_center), 
                    m_params_BosonStar(a_params_BosonStar), 
                    m_params_BosonStar2(a_params_BosonStar2)
    {   
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        weightfunc_t<data_t> out;

        Coordinates<double> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);

        Coordinates<double> coords_grid(current_cell, m_dx,
        m_center);

        double rapidity = m_params_BosonStar.BS_rapidity;
        double rapidity2 = m_params_BosonStar2.BS_rapidity;
        double separation = m_params_BosonStar.BS_separation;
        double impact_parameter = m_params_BosonStar.BS_impact_parameter;

        double c_1 = cosh(rapidity);
        double x1 = (coords.x-separation/2.) * c_1;
        double y1 = coords.y+impact_parameter/2.;

        double argument1 = (1/separation) * (sqrt(pow(coords_grid.x-x1, 2)+pow(coords_grid.y-y1,2)));

        double c_2 = cosh(-rapidity2);
        double x2 = (coords.x+separation/2.)*c_2;
        double y2 = coords.y-impact_parameter/2.;

        double argument2 = (1/separation) * (sqrt(pow(coords_grid.x-x2, 2)+pow(coords_grid.y-y2,2)));

        double weight_func1 = 42.0;
        double weight_func2 = 42.0;

        WeightFunction weightfunction;

        if ((coords_grid.z = m_center[2]))
        {
            weight_func1 = weightfunction.compute_weight(argument1); // bump at object 1
            weight_func2 = weightfunction.compute_weight(argument2); //bump at object 2
        }
        else
        {
            weight_func1  = 0.0;
            weight_func2 = 0.0;
        }

    out.weight1 = weight_func1;
    out.weight2 = weight_func2;

    current_cell.store_vars(out.weight1, c_weight1);
    current_cell.store_vars(out.weight2, c_weight2);
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */