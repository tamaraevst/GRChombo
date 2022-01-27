#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

class WeightFunction
{
    public:

    WeightFunction(){}

    double compute_weight(double scaledr) const
    {   
        double weightfunc;
      
        if (scaledr <= 1.0)
        { 
            return weightfunc = 140.0 * ( (1.0/4.0)*pow((1-scaledr),4) - (3.0/5.0)*pow((1-scaledr),5) + (1.0/2.0)* pow((1-scaledr),6) - (1.0)/(7.0)*pow((1-scaledr),7));
            //return weightfunc = 6.0 * ((1.0/2.0) * pow((1-scaledr), 2) - (1.0/3.0) * pow((1-scaledr), 3));
	//return weightfunc = 4808643120 * ((1.0/16.0) * pow(1-scaledr, 16) - (15.0/17.0) * pow(1-scaledr, 17) + (35.0/6.0) * pow(1-scaledr, 18) - (455.0/19.0) * pow(1-scaledr, 19) + (273.0/4.0) * pow(1-scaledr, 20) - 143.0 * pow(1-scaledr, 21) + (455.0/2.0) * pow(1-scaledr, 22) - (6435.0/23.0) * pow(1-scaledr, 23) + (2145.0/8.0) * pow(1-scaledr, 24) - (1001.0/5.0) * pow(1-scaledr, 25) + (231.0/2.0) * pow(1-scaledr, 26) - (455.0/9.0) * pow(1-scaledr, 27) + (65.0/4.0) * pow(1-scaledr, 28) - (105.0/29.0) * pow(1-scaledr, 29) + (1.0/2.0) * pow(1-scaledr, 30) - (1.0/31.0) * pow(1-scaledr, 31));  
	}

        else
        {
            return weightfunc = 0.0;
        }

    }

    double stretching_factor(double coord_x, double coord_y, double alpha_zero)
    {
        double sin_half_angle_squared = (1.0 - (coord_x) / sqrt(pow(coord_x, 2) + pow(coord_y, 2))) / (2.0);
        
        double stretch = alpha_zero + (1.0 - alpha_zero) * sin_half_angle_squared;

        return stretch;
    }
    
};

#endif /* WEIGHTFUNCTION_HPP_ */
