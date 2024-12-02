/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SINGLEBOSONSTAR_HPP_)
#error "This file should only be included through SingleBosonStar.hpp"
#endif

#ifndef SINGLEBOSONSTAR_IMPL_HPP_
#define SINGLEBOSONSTAR_IMPL_HPP_

#include <iostream>
#include <vector>
#include <fstream> // Needed for std::ifstream

SingleBosonStar::SingleBosonStar() {}

// Read below 
void SingleBosonStar::main(std::string m_base_path)
{
    std::vector<double> temp_data_x1;  // array for theta & array for radius 
    double value_x1;
    // std::string base_path = "/Users/macbookpro/Desktop/PHD/Projects/InstabilityRing/";
    std::string filename_radius = "radius.dat";
    std::string filename_amp = "urA.dat";
    std::string filename_f = "f.dat";
    std::string filename_lapse = "urPhi.dat";

    std::string full_path_radius = m_base_path + filename_radius;
    std::string full_path_amp = m_base_path + filename_amp;
    std::string full_path_f = m_base_path + filename_f;
    std::string full_path_lapse = m_base_path + filename_lapse;

    ////////////
    // Radius //
    ////////////

    std::ifstream file_x1(full_path_radius);
    if (file_x1.is_open()) {
        // Read each value from the file
        while (file_x1 >> value_x1) {
            temp_data_x1.push_back(value_x1);  // Append to temporary vector
        }
        file_x1.close();
    } else {
        std::cerr << "Unable to open file radius.dat" << std::endl;
    }
    m = temp_data_x1.size();
    x1.resize(m);
    for (int i = 0; i < m; ++i) {
        x1[i] = temp_data_x1[i];  // Copy from std::vector to VecDoub
    }

    // Resize all matrix variables
    amp.resize(m);
    f.resize(m);     
    Phi.resize(m);    

    ///////////////
    // Amplitude //
    ///////////////

    std::ifstream file_A(full_path_amp);
    if (file_A.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_A >> amp[i];
            }
        file_A.close();
    } else {
            std::cerr << "Unable to open file urA.dat" << std::endl;
        }

    /////////////////////
    // f conformal var //
    /////////////////////

    std::ifstream file_f(full_path_f);
    if (file_f.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_f >> f[i];
        }
        file_f.close();
        } else {
            std::cerr << "Unable to open file f.dat" << std::endl;
        }

    ////////////
    // lapse //
    ///////////

    std::ifstream file_l(full_path_lapse);
    if (file_l.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_l >> Phi[i];
        }
        file_l.close();
        } else {
            std::cerr << "Unable to open file urPhi.dat" << std::endl;
        }

}

// double SingleBosonStar::get_amp_interp(double r, int m) const
// {
//     double interpolated_value; 
//     if (x1.size()==0) {
//     MayDay::Error("Error: x1 is empty. Make sure to initialize it before calling get_amp_interp().");
//     }

//     Poly_interp cubic_inerpolation(x1,amp,4); 
//     interpolated_value = cubic_inerpolation.interp(r);
    
//     return interpolated_value;
// }

// double SingleBosonStar::get_f_interp(const double r, int m) const
// {
//     double interpolated_value; 

//     Poly_interp cubic_inerpolation(x1,f,4); 
//     interpolated_value = cubic_inerpolation.interp(r);
    
//     return interpolated_value;
// }

// double SingleBosonStar::get_lapse_interp(const double r, int m) const
// {
//     double interpolated_value; 

//     Poly_interp cubic_inerpolation(x1,Phi,4); 
//     interpolated_value = cubic_inerpolation.interp(r);
    
//     return interpolated_value;
// }

double SingleBosonStar::get_amp_interp(const double r) const
{   
    double dx = 0.00099974709747009;
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? amp[1] : amp[iter - 1]); // conditionl/ternary imposing
                                               // zero gradeint at r=0
    f2 = amp[iter];
    f3 = amp[iter + 1];
    f4 = amp[iter + 2];

    if (iter > 1000001 - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double SingleBosonStar::get_f_interp(const double r) const
{   
    double dx = 0.00099974709747009;
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? f[1] : f[iter - 1]); // conditionl/ternary imposing
                                               // zero gradeint at r=0
    f2 = f[iter];
    f3 = f[iter + 1];
    f4 = f[iter + 2];

    if (iter > 1000001 - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

double SingleBosonStar::get_lapse_interp(const double r) const
{   
    double dx = 0.00099974709747009;
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? Phi[1] : Phi[iter - 1]); // conditionl/ternary imposing
                                               // zero gradeint at r=0
    f2 = Phi[iter];
    f3 = Phi[iter + 1];
    f4 = Phi[iter + 2];

    if (iter > 1000001 - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}

#endif /* SINGLEBOSONSTAR_IMPL_HPP_ */