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
        }

        else
        {
            return weightfunc = 0.0;
        }

    }

    
};

#endif /* WEIGHTFUNCTION_HPP_ */