#include "FArrayBox.H"
#include <iostream>
#include "UserVariables.hpp"
#include "SetValue.hpp"
#include "tensor.hpp"
#include "GRutils.hpp"

#include "BoxLoops.hpp"
#include "SetValue.hpp"



int main()
{
   int passed = 1;

    const int N_GRID = 8;
    Box box(IntVect(0,0,0), IntVect(N_GRID-1,N_GRID-1,N_GRID-1));
    FArrayBox in_fab(box, c_NUM);

    double value = 42.;
    BoxLoops::loop(SetValue(42.),in_fab,in_fab);

    for (int iz = 0; iz < N_GRID; ++iz)
    {
        for (int iy = 0; iy < N_GRID; ++iy)
        {
            for (int ix = 0; ix < N_GRID; ++ix)
            {
                const IntVect iv(ix,iy,iz);
                if (in_fab(iv, c_chi) != value) 
                {
                    passed = -1;
                    pout () << iv << std::endl;
                }
            }
        }
    }


   if (passed == 1) std::cout << "SetValue test passed" << std::endl;
   else std::cout << "SetValue test NOT passed" << std::endl;

   return passed;
}
