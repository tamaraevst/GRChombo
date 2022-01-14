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

        double argument1 = (1.0/separation) * (sqrt(pow(coords.x-separation*q/(q+1), 2)+pow(coords.y,2)+pow(coords.z, 2)));

        double argument2 = (1.0/ separation) * (sqrt(pow(coords.x+separation/(q+1), 2)+pow(coords.y,2)+pow(coords.z, 2)));
	
        double weight_func1 = 42.0;
        double weight_func2 = 42.0;

        WeightFunction weightfunction;

        weight_func1 = weightfunction.compute_weight(argument1); // bump at object 1
        weight_func2 = weightfunction.compute_weight(argument2); //bump at object 2
       
	
    	out.weight1 = weight_func1;
    	out.weight2 = weight_func2;

    	current_cell.store_vars(out.weight1, c_weight1);
    	current_cell.store_vars(out.weight2, c_weight2);
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */
