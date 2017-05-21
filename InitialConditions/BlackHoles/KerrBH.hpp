//Last update K Clough 21.05.2017

#ifndef KERRBH_HPP_
#define KERRBH_HPP_

#include "simd.hpp"
#include "VarsBase.hpp"
#include "tensor.hpp"
#include "Coordinates.hpp"
#include "CCZ4.hpp"
#include "Cell.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

//! Class which computes the Kerr initial conditions per arXiv 1401.1548
class KerrBH
{
//Use the variable definition in CCZ4
template<class data_t>
using Vars=CCZ4::Vars<data_t>;

public:
    //! Stuct for the params of the Kerr BH
    struct params_t
    {
        double mass;//!<< The mass of the Kerr BH
        std::vector<double> center;//!< The center of the Kerr BH
        double spin;//!< The spin param a = J/M, so 0 < a < 1
    };

protected:
    double m_dx;
    params_t m_params;

public:
    KerrBH(params_t a_params, double a_dx) : m_dx (a_dx), m_params (a_params)

    {
        // check this spin param is sensible
        if((m_params.spin > 1.0) || (m_params.spin < 0.0))
        {
            MayDay::Error("The spin parameter must be in the range 0 < a < 1.0");
        }
    }

    template <class data_t>
    void compute(Cell current_cell);

protected:
    //! Function which computes the components of the metric in spherical coords
    template <class data_t>
    void compute_kerr(tensor<2,data_t> &spherical_g,//!<< The spatial metric in spherical coords
                      tensor<2,data_t> &spherical_K,//!<< The extrinsic curvature in spherical coords
                      tensor<1,data_t> &spherical_shift,//!<< The spherical components of the shift
                      data_t &kerr_lapse,//!<< The lapse for the kerr solution
                      const Coordinates<data_t> coords);//!<< Coords of current cell

    //! Function which computes the cartesian components of the CCZ4 vars from spherical K_ij, g_ij and beta^i
    template <class data_t>
    void spherical_to_cartesian(Vars<data_t> &vars,
                                tensor<2,data_t> &spherical_g,//!<< The spatial metric in spherical coords
                                tensor<2,data_t> &spherical_K,//!<< The extrinsic curvature in spherical coords
                                tensor<1,data_t> &spherical_shift, //!<< The components of the kerr shift in spherical coords
                                const Coordinates<data_t> coords);//!<< Coords of current cell

    //! Function to get coordinates in cartesian space, and radius
    template <class data_t>
    void get_position(Coordinates<data_t> coords, data_t &x, double &y, double &z, 
                          data_t &r, data_t &r2, data_t &rho, data_t &rho2);
};

#include "KerrBH.impl.hpp"

#endif /* KERRBH_HPP_ */
