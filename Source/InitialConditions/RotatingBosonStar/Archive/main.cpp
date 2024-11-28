#include <iostream>
#include "RotatingBosonStarSolution.hpp"
#include "RotatingBosonStarSolution.impl.hpp"
#include "cmath"
#include <vector>

int main()
{   
	RotatingBosonStarSolution a_boson_star;

    a_boson_star.read();
    double ff = a_boson_star.get_BSfrequency();
    std::cout << "BS frequency is ..." << ff << std::endl;

    double r = 0.501002;
    double theta = 3.07845;

    double dthomega_val = a_boson_star.get_dthomega_interp(r, theta);
    std::cout << "Dthomega value is ..." << dthomega_val << std::endl;

    double A_val = a_boson_star.get_amp_interp(r, theta);
    std::cout << "A value is ..." << A_val << std::endl;

    double f_val = a_boson_star.get_f_interp(r, theta);
    std::cout << "f value is ..." << f_val << std::endl;

	return 0;
}

