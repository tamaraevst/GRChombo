#ifdef _OPENMP
#include <omp.h>
#endif

#include "FArrayBox.H"
#include "UserVariables.hpp"
#include <sys/time.h>
#include "ScalarField.hpp"
#include "HarmonicTest.hpp"
#include "Cell.hpp"
#include "BoxLoops.hpp"
#include "DebuggingTools.hpp"
#include "ComputeClassPack.hpp"
#include "SetValue.hpp"

int main()
{
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    const int N_GRID = 64;
    Box box(IntVect(0,0,0), IntVect(N_GRID-1,N_GRID-1,N_GRID-1));
    Box ghosted_box(IntVect(-3,-3,-3), IntVect(N_GRID+2,N_GRID+2,N_GRID+2));
    FArrayBox in_fab(ghosted_box, c_NUM);
    in_fab.setVal(0.0);
//    BoxLoops::loop(make_compute_pack(SetValue(0.0)), in_fab, in_fab);
    FArrayBox out_fab(box, c_NUM);
    out_fab.setVal(0.0);
//    BoxLoops::loop(make_compute_pack(SetValue(0.0)), out_fab, out_fab);
    double length = 64.0;

    const double dx = length / (N_GRID);
    const double center = length/2.0;

    for (int iz = -3; iz < N_GRID+3; ++iz)
    {
        const double z = (iz+0.5)*dx - center;
        for (int iy = -3; iy < N_GRID+3; ++iy)
        {
            const double y = (iy+0.5)*dx - center;
            for (int ix = -3; ix < N_GRID+3; ++ix)
            {
                const double x = (ix+0.5)*dx - center;
                double r = sqrt(x*x + y*y + z*z);
                double rho = sqrt(x*x + y*y);
                const IntVect iv(ix,iy,iz);
                if (r < 1e-6)
                {
                    r = 1e-6;
                }
                if (rho < 1e-6)
                {
                    rho = 1e-6;
                }

                // here testing the es = -1, el = 2, em = -1 case
                // and also the calculation of r in coords
                double harmonic;
                harmonic = sqrt(5.0/16.0/M_PI)*x*(2*z*z - z*r - r*r)/rho/r/r;
                in_fab(iv, c_phi) = harmonic/r;
            }
        }
    }

    std::vector<double> center_vector = {center, center, center};

    //Test the spherical harmonics across grid
    BoxLoops::loop(HarmonicTest(center_vector, dx), in_fab, out_fab); // disable_simd());
    out_fab -= in_fab;

    for (int i = 0; i < c_NUM; ++i)
    {
        double max_err = out_fab.norm(0, i, 1);
        double max_act = in_fab.norm(0,i,1);
        if (max_err/max_act > 1e-10)
        {
            std::cout << "COMPONENT " << UserVariables::variable_names[i] << " DOES NOT AGREE: MAX ERROR = " << out_fab.norm(0, i, 1) << std::endl;
            std::cout << "COMPONENT " << UserVariables::variable_names[i] << " DOES NOT AGREE: MAX Actual Value = " << max_act << std::endl;
        }
    }
}
