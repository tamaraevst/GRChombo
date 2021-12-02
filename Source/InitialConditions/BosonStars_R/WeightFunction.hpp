#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

class WeightFunction
{
    public:

    double weightfunction(double scaledr) const
    {
        double weightfunc = 140 * ( (1.0/4.0)*pow((1-scaledr),4) - (3.0/5.0)*pow((1-scaledr),5) + (1.0/2.0)* pow((1-scaledr),6) - (1.0)/(7.0)*pow((1-scaledr),7));

        if (scaledr <= 1)
        {
            return weightfunc;
        }

        if (scaledr > 1)
        {
            return 0;
        }
    }

};

#endif /* WEIGHTFUNCTION_HPP_ */