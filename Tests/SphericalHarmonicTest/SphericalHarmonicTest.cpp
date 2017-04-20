#ifdef _OPENMP
#include <omp.h>
#endif

#include "FArrayBox.H"
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include "FABDriver.hpp"
#include "ScalarField.hpp"
#include "HarmonicTest.hpp"
#include "UserVariables.hpp"
#include "DebuggingTools.hpp"

#define CHF_FRAn(a, n, c) \
    a.dataPtr(n), \
    D_DECL6(&a.loVect()[0],  \
            &a.loVect()[1],  \
            &a.loVect()[2],  \
            &a.loVect()[3],  \
            &a.loVect()[4],  \
            &a.loVect()[5]), \
    D_DECL6(&a.hiVect()[0],  \
            &a.hiVect()[1],  \
            &a.hiVect()[2],  \
            &a.hiVect()[3],  \
            &a.hiVect()[4],  \
            &a.hiVect()[5]), \
    &c

#define CHF_CONST_FRAn(a, n, c) CHF_FRAn(a ,n, c)

int main()
{
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif

    const int N_GRID = 64;
    Box box(IntVect(0,0,0), IntVect(N_GRID-1,N_GRID-1,N_GRID-1));
    Box ghosted_box(IntVect(-3,-3,-3), IntVect(N_GRID+2,N_GRID+2,N_GRID+2));
    FArrayBox in_fab(ghosted_box, c_NUM);
    in_fab.setVal(0.);
    FArrayBox out_fab(box, c_NUM);
    out_fab.setVal(0.);
    double length = 64.0;

    const double dx = length / (N_GRID);
    const double center = length/2.0;

    for (int zz = -3; zz < N_GRID+3; ++zz)
    {
        const double z = (zz+0.5)*dx - center;
        for (int yy = -3; yy < N_GRID+3; ++yy)
        {
            const double y = (yy+0.5)*dx - center;
            for (int xx = -3; xx < N_GRID+3; ++xx)
            {
                const double x = (xx+0.5)*dx - center;
                double r = sqrt(x*x + y*y + z*z);
                double rho = sqrt(x*x + y*y);
                const IntVect iv(xx,yy,zz);
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
//                if ((xx==30) && (yy==30) && (zz==30))
//                {
//                    DEBUG_OUT(harmonic);
//                }
            }
        }
    }

    //Fill in the Parameters
//    std::vector<double> center_vector;
    IntVect center_vector(center, center, center);

    //Test the spherical harmonics across grid
    FABDriver<HarmonicTest>(center_vector, dx).execute(in_fab, out_fab); // disable_simd());
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
