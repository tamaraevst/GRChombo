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
        bool do_stretch = m_params_BosonStar.do_stretch;
        int n_weight = m_params_BosonStar.n_power;
        int id_choice = m_params_BosonStar.id_choice;
        double radius_width = m_params_BosonStar.radius_width;

        WeightFunction weightfunction;

        if (id_choice == 2)
        {
            double factor1, factor2;

            if (do_stretch)
            {
                factor1 = weightfunction.stretching_factor((coords.x-separation/(q+1))*cosh(rapidity), coords.y, alpha);
                factor2 = weightfunction.stretching_factor2((coords.x+q*separation/(q+1))*cosh(-rapidity2), coords.y, alpha);
            }
            else
            {
                factor1 = 1.0;
                factor2 = 1.0;
            }

            double argument1 = (factor1/separation) * (sqrt(pow((coords.x-separation/(q+1))*cosh(rapidity), 2)+pow(coords.y,2)+pow(coords.z, 2)));
            double argument2 = (factor2/ separation) * (sqrt(pow((coords.x+q*separation/(q+1))*cosh(-rapidity2), 2)+pow(coords.y,2)+pow(coords.z, 2)));
	
            double weight_func1 = 42.0;
            double weight_func2 = 42.0;
        
            weight_func1 = weightfunction.compute_weight(argument1, n_weight); // bump at object 1
            weight_func2 = weightfunction.compute_weight(argument2, n_weight); //bump at object 2
       
    	    out.weight1 = weight_func1;
    	    out.weight2 = weight_func2;

    	    current_cell.store_vars(out.weight1, c_weight1);
    	    current_cell.store_vars(out.weight2, c_weight2);
        }

        if (id_choice == 3)
        {
            double profile_func1 = weightfunction.profile_chi((coords.x-separation/(q+1))*cosh(rapidity), coords.y, coords.z, radius_width);
            double profile_func2 = weightfunction.profile_chi((coords.x+q*separation/(q+1))*cosh(-rapidity2), coords.y, coords.z, radius_width);

            current_cell.store_vars(profile_func1, c_profile1);
            current_cell.store_vars(profile_func2, c_profile2);
        }
    
    }

    
};

#endif /* COMPUTEWEIGHTFUNCTION_HPP_ */
