#include "FArrayBox.H"
#include <iostream>
#include "PositiveChiAndAlpha.hpp"
#include "tensor.hpp"
#include "GRutils.H"

#include "FABDriver.hpp"
#include "PositiveChiAndAlpha.hpp"



int main()
{
   int passed = 1;

    const int N_GRID = 8;
    Box box(IntVect(0,0,0), IntVect(N_GRID-1,N_GRID-1,N_GRID-1));
    FArrayBox in_fab(box, c_NUM);

    for (int iz = 0; iz < N_GRID; ++iz)
    {
        for (int iy = 0; iy < N_GRID; ++iy)
        {
            for (int ix = 0; ix < N_GRID; ++ix)
            {
                const IntVect iv(ix,iy,iz);
                double value;
                if (ix < N_GRID/2) value = 1;
                else value = 1e-10;

                in_fab(iv, c_chi) = value;
                in_fab(iv, c_lapse) = value;
            }
        }
    }

    FABDriver<PositiveChiAndAlpha>().execute(in_fab,in_fab);

    for (int iz = 0; iz < N_GRID; ++iz)
    {
        for (int iy = 0; iy < N_GRID; ++iy)
        {
            for (int ix = 0; ix < N_GRID; ++ix)
            {
                const IntVect iv(ix,iy,iz);
                double value;
                if (ix < N_GRID/2) value = 1; //PositiveChiAndAlpha should leave this untouched
                else value = 1e-4;  //PositiveChiAndAlpha should change 1e-10 to 1e-4

                if ( (in_fab(iv, c_chi) != value)
                     || (in_fab(iv, c_lapse) != value) ) passed = -1;
            }
        }
    }


   if (passed == 1) std::cout << "PositiveChiAndAlpha test passed" << std::endl;
   else std::cout << "PositiveChiAndAlpha test NOT passed" << std::endl;

   return passed;
}
