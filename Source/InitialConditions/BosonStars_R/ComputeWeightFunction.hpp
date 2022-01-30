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
        double alpha = m_params_BosonStar.alpha_stretch;

        double c_ = cosh(rapidity);
        double x = (coords.x-separation/(q+1))*c_;
        double y = coords.y+impact_parameter/2.;
        double z = coords.z; //set /tilde{t} to zero

        double c_2 = cosh(-rapidity2);
        double x2 = (coords.x+q*separation/(q+1))*c_2;
        double y2 = coords.y-impact_parameter/2.;

        WeightFunction weightfunction;

        double factor1 = weightfunction.stretching_factor(x, y, alpha);
        double argument1 = (factor1/separation) * (sqrt(pow(x, 2)+pow(y,2)+pow(z, 2)));

        double factor2 = weightfunction.stretching_factor2(x2, y2, alpha);
        double argument2 = (factor2/ separation) * (sqrt(pow(x2, 2)+pow(y2,2)+pow(z, 2)));
	
        double weight_func1 = 42.0;
        double weight_func2 = 42.0;
        
        weight_func1 = weightfunction.compute_weight(argument1); // bump at object 1
        weight_func2 = weightfunction.compute_weight(argument2); //bump at object 2
       
    	out.weight1 = weight_func1;
    	out.weight2 = weight_func2;

    	current_cell.store_vars(out.weight1, c_weight1);
    	current_cell.store_vars(out.weight2, c_weight2);
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */
