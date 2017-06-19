#ifndef COMBINATORICS_HPP_
#define COMBINATORICS_HPP_

namespace Combinatorics
{
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
        CH_assert( (r>=0) && (n>=r) ); //sense check values

        double out = factorial(n) / (factorial(r) * factorial(n-r));
        return out;
    }
}

#endif /* COMBINATORICS_HPP_ */
