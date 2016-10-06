#ifndef BUBBLESF_HPP_
#define BUBBLESF_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

class BubbleSF
{
protected:
    const FABDriverBase& m_driver;
    double m_dx;

public:
    struct params_t
    {
        Real amplitudeSF;
        std::vector<Real> centerSF;
        Real widthSF;
    };

    const params_t m_params;

    BubbleSF(const FABDriverBase& a_driver, params_t a_params,
             double a_dx) :
        m_driver (a_driver), m_dx (a_dx), m_params (a_params) {}

    //Not currently vectorised (it is only done once so it's hardly worth adding all the special functions to simd)
    void compute(int ix, int iy, int iz);

protected:
    double compute_phi(double x, double y, double z);

};

void BubbleSF::compute(int ix, int iy, int iz)
{
    CCZ4SFMatter::vars_t<double> vars;
    vars.assign(0.); //Set only the non-zero components explicitly below
    Coordinates<double> coords(ix,iy,iz,m_dx);
    double x = coords.x; //TODO: change functions to accept coords rather than x,y,z
    double y = coords.y;
    double z = coords.z;

    vars.phi = compute_phi(x,y,z);
    vars.PiM = 0;

    vars.lapse = 1;

    vars.chi = 1;
    //Conformal metric is flat
    FOR1(i) vars.h[i][i] = 1.;

    m_driver.store_vars(vars);

}

double BubbleSF::compute_phi(double x, double y, double z)
{
	  Real rr2 = pow(x - m_params.centerSF[0],2) + pow(y - m_params.centerSF[1],2) + pow(z - m_params.centerSF[2],2);

  	if (std::fabs(rr2) < 1e-6)
  	{
	 	  rr2 = 1e-6;
  	}

    const double phiout = m_params.amplitudeSF*exp(-rr2/m_params.widthSF);
    return phiout;
}

#endif /* BUBBLESF_HPP_ */
