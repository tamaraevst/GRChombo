#ifndef SPHERICALHARMONICS_HPP_
#define SPHERICALHARMONICS_HPP_

#include "tensor.hpp"
#include "Coordinates.hpp"

// Functions for the spin weighted spherical harmonics
// The harmonics here are taken from recurrance relation for Wigner small d-matrix elements
// https://en.wikipedia.org/wiki/Wigner_D-matrix (NB THE WIKIPEDIA FORMULA IS NOT QUITE CORRECT)
// https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics
namespace SphericalHarmonics
{
    template <class data_t>
    struct Y_lm_t {
        data_t Real;
        data_t Im;
    };

    //Calculate factorials
    double factorial(int n)
    {
      double out = 1.0;
      for (int i = 1; i <= n; i++)
      {
          out *= i;
      }
      return out;
    }

    //Calculate combinatorics
    double n_choose_r(int n, int r)
    {
        //sense check values
        if ((n < 0) || (r > n) || (r < 0))
        {
            //Invalid values for combinatorics
            MayDay::Error("Nonsensical choice of n and r in combinatorics, in utils/SphericalHarmonics.hpp");
        }

        double out = factorial(n) / (factorial(r) * factorial(n-r));

        return out;
    }

    // Calculates the spin weight es, el, em spherical harmonic
    template <class data_t>
    Y_lm_t<data_t> spin_Y_lm(const data_t x,
                              const double y,
                              const double z,
                              const data_t r,
                              const int es,
                              const int el,
                              const int em)
    {
        if ((el < abs(em)) || (el < 0))
        {
            //Invalid values for el and em
            MayDay::Error("Nonsensical choice of el and em, in utils/SphericalHarmonics.hpp");
        }

        Y_lm_t<data_t> Y_lm; 

        data_t theta = acos(z/r);
        data_t phi = atan2(y,x);

        double coefficient = pow(-1.0, em) * sqrt((2.0*el + 1.0)/(4.0*M_PI));
        coefficient *= sqrt(factorial(el+em)*factorial(el-em)/factorial(el+es)/factorial(el-es));

        data_t sum = 0.0;
        int lower_limit = em-es > 0 ? em-es : 0;
        int upper_limit = el+em < el-es ? el+em : el-es;

        for(int i = lower_limit; i <= upper_limit; i++)
        {
            double temp = n_choose_r(el - es, i) * n_choose_r(el + es, i + es - em);
            sum += temp * pow(-1.0, el - i - es)*pow(cos(theta/2.0), 2*i + es - em)*pow(sin(theta/2.0), 2*(el-i) + em - es);
        }

        Y_lm.Real = coefficient*sum*cos(em*phi);
        Y_lm.Im   = coefficient*sum*sin(em*phi);

        return Y_lm;
    }

}

#endif /* SPHERICALHARMONICS_HPP_ */
