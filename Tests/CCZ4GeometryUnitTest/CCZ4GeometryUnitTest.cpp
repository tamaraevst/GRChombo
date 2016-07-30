#include <iostream>
#include "CCZ4Geometry.hpp"
#include "tensor.hpp"
#include "GRutils.H"

template <class data_t>
struct vars_t
{
   data_t chi;
   tensor<2, data_t> h;
   tensor<1, data_t> Gamma;
};


int main()
{
   int passed = 1;

   vars_t<double> vars;
   vars_t< tensor<1,double> > d1;
   vars_t< tensor<2,double> > d2;
   tensor<1, double>  Z_over_chi;

#include "values1.hpp" //Including the auto generated file with values

   auto h_UU = TensorAlgebra::compute_inverse(vars.h);

   auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

   auto ricciZ = CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

   //Compare
   FOR2(i,j)
   {
      double diff = h_UU[i][j] - h_UU_known[i][j];
      if (diff > 1e-14)
      {
         std::cout << "h_UU wrong in component [" << i << "]["<< j<< "]"<< std::endl;
         passed = -1;
      }
   }

   FOR3(i,j,k)
   {
      double diff = chris.ULL[i][j][k] - chris_known[i][j][k];
      if (diff > 1e-14)
      {
         std::cout << "chris wrong in component [" << i << "]["<< j<< "][" << k << "]" << std::endl;
         std::cout << "value: " << chris.ULL[i][j][k] << std::endl;
         std::cout << "correct value: " << chris_known[i][j][k] << std::endl;
         passed = -1;
      }
   }

   FOR1(i)
   {
      double diff = chris.contracted[i] - chris_contracted_known[i];
      if (diff > 1e-14)
      {
         std::cout << "chris contracted wrong in component [" << i << "]" << std::endl;
         std::cout << "value: " << chris.contracted[i] << std::endl;
         std::cout << "correct value: " << chris_contracted_known[i] << std::endl;
         passed = -1;
      }
   }

   FOR2(i,j)
   {
      double diff = ricciZ.LL[i][j] - ricciZ_known[i][j];
      if (diff > 1e-14)
      {
         std::cout << "ricciZ contracted wrong in component [" << i << "][" << j << "]" << std::endl;
         std::cout << "value: " << ricciZ.LL[i][j] << std::endl;
         std::cout << "correct value: " << ricciZ_known[i][j] << std::endl;
         passed = -1;
      }
   }

   double diff = ricciZ.scalar - ricciZ_scalar_known;
   if (diff > 1e-14)
   {
      std::cout << "ricci scalar wrong" << std::endl;
      std::cout << "value: " << ricciZ.scalar << std::endl;
      std::cout << "correct value: " << ricciZ_scalar_known << std::endl;
      passed = -1;
   }


   if (passed == 1) std::cout << "CCZ4Geometry test passed" << std::endl;
   else std::cout << "CCZ4Geometry test NOT passed" << std::endl;

   return passed;
}
