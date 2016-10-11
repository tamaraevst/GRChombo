#include <iostream>
#include "simd.hpp"

int main()
{
   int error = 0;

   int simd_length = simd_traits<double>::simd_len;

   //Test the log
   double* xLog = new double[simd_length];
   for (int i = 0; i < simd_length; ++i) xLog[i] = exp(i+1);
   auto simd_in = simd<double>::load(xLog);
   auto simd_out = simd_log(simd_in);
   simd<double>::store(xLog, simd_out);
   for (int i = 0; i < simd_length; ++i) if (xLog[i] != i+1) error = -1;
   delete[] xLog;

   //Test the sqrt
   double* xSqrt = new double[simd_length];
   for (int i = 0; i < simd_length; ++i) xSqrt[i] = pow(i+1,2);
   simd_in = simd<double>::load(xSqrt);
   simd_out = simd_sqrt(simd_in);
   simd<double>::store(xSqrt, simd_out);
   for (int i = 0; i < simd_length; ++i) if (xSqrt[i] != i+1) error = -2;
   delete[] xSqrt;

   if (error == 0) std::cout << "Simd functions unit test passed" << std::endl;
   else std::cout << "Simd functions unit test NOT passed" << std::endl;

   return error;
}
