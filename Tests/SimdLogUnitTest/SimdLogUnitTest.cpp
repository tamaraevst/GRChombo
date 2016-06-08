#include "../../simd/simd.hpp"
#include <iostream>

int main()
{
   double x[4] = { 7.389, 2.7182818, 7.389, 2.7182818 };
   auto simd_in = simd<double>::load(x);
   auto simd_out = simd_log(simd_in);
   simd<double>::store(x, simd_out);
   for (int i = 0; i < 4; ++i)
      std::cout << x[i] << std::endl;
}
