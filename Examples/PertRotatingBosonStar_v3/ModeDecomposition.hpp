/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 #ifndef MODEDECOMPOSITION_HPP
 #define MODEDECOMPOSITION_HPP
 
 #include "ComplexScalarField.hpp"
 #include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
 #include "Cell.hpp"
 #include "Coordinates.hpp"
 #include "UserVariables.hpp"
 #include "FourthOrderDerivatives.hpp"
 #include "simd.hpp"
 
//  template <class deriv_t = FourthOrderDerivatives>
 class ModeDecomposition
 {
   protected:
    const FourthOrderDerivatives m_deriv;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center

    template <class data_t> using ADMVars
                                = ADMConformalVars::VarsNoGauge<data_t>;
    template <class data_t> using MatterVars
                                = ComplexScalarField<>::Vars<data_t>;
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

   public:
 
    //! Constructor
    ModeDecomposition(const double a_dx, const std::array<double, CH_SPACEDIM> a_center) : m_dx(a_dx), m_center(a_center), m_deriv(a_dx) {}
    
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto vars = current_cell.template load_vars<Vars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto advec_csf = m_deriv.template advection<MatterVars>(current_cell, vars.shift);

        Coordinates<double> coords(current_cell, m_dx, m_center);
        double y = coords.y;
        double x = coords.x;
        double phi = atan2(y,x);
                      
        //Calculate time derivative phi^2
        data_t phi_Re_t = advec_csf.phi_Re - vars.lapse * matter_vars.Pi_Re;
        data_t phi_Im_t = advec_csf.phi_Im - vars.lapse * matter_vars.Pi_Im;
        data_t dt_A_sq = abs(2 * phi_Re_t * matter_vars.phi_Re + 2 * phi_Im_t * matter_vars.phi_Im);

        data_t mode_dtA_01_re = dt_A_sq * cos(1*phi);
        data_t mode_dtA_01_im = dt_A_sq * sin(1*phi);
        data_t mode_dtA_02_re = dt_A_sq * cos(2*phi);
        data_t mode_dtA_02_im = dt_A_sq * sin(2*phi);
        data_t mode_dtA_03_re = dt_A_sq * cos(3*phi);
        data_t mode_dtA_03_im = dt_A_sq * sin(3*phi);
        data_t mode_dtA_04_re = dt_A_sq * cos(4*phi);
        data_t mode_dtA_04_im = dt_A_sq * sin(4*phi);
        data_t mode_dtA_05_re = dt_A_sq * cos(5*phi);
        data_t mode_dtA_05_im = dt_A_sq * sin(5*phi);
        data_t mode_dtA_06_re = dt_A_sq * cos(6*phi);
        data_t mode_dtA_06_im = dt_A_sq * sin(6*phi);
        data_t mode_dtA_07_re = dt_A_sq * cos(7*phi);
        data_t mode_dtA_07_im = dt_A_sq * sin(7*phi);
        data_t mode_dtA_08_re = dt_A_sq * cos(8*phi);
        data_t mode_dtA_08_im = dt_A_sq * sin(8*phi);
        data_t mode_dtA_09_re = dt_A_sq * cos(9*phi);
        data_t mode_dtA_09_im = dt_A_sq * sin(9*phi);
        data_t mode_dtA_10_re = dt_A_sq * cos(10*phi);
        data_t mode_dtA_10_im = dt_A_sq * sin(10*phi);
        data_t mode_dtA_11_re = dt_A_sq * cos(11*phi);
        data_t mode_dtA_11_im = dt_A_sq * sin(11*phi);
        data_t mode_dtA_12_re = dt_A_sq * cos(12*phi);
        data_t mode_dtA_12_im = dt_A_sq * sin(12*phi); 
        data_t mode_dtA_13_re = dt_A_sq * cos(13*phi);
        data_t mode_dtA_13_im = dt_A_sq * sin(13*phi);
        data_t mode_dtA_14_re = dt_A_sq * cos(14*phi);
        data_t mode_dtA_14_im = dt_A_sq * sin(14*phi);
        data_t mode_dtA_15_re = dt_A_sq * cos(15*phi);
        data_t mode_dtA_15_im = dt_A_sq * sin(15*phi);
        data_t mode_dtA_16_re = dt_A_sq * cos(16*phi);
        data_t mode_dtA_16_im = dt_A_sq * sin(16*phi);
        data_t mode_dtA_17_re = dt_A_sq * cos(17*phi);
        data_t mode_dtA_17_im = dt_A_sq * sin(17*phi);
        data_t mode_dtA_18_re = dt_A_sq * cos(18*phi);
        data_t mode_dtA_18_im = dt_A_sq * sin(18*phi);
        data_t mode_dtA_19_re = dt_A_sq * cos(19*phi);
        data_t mode_dtA_19_im = dt_A_sq * sin(19*phi);
        data_t mode_dtA_20_re = dt_A_sq * cos(20*phi);
        data_t mode_dtA_20_im = dt_A_sq * sin(20*phi);
                    
        current_cell.store_vars(mode_dtA_01_re, c_mode_dtA_01_re);
        current_cell.store_vars(mode_dtA_01_im, c_mode_dtA_01_im);
        current_cell.store_vars(mode_dtA_02_re, c_mode_dtA_02_re);
        current_cell.store_vars(mode_dtA_02_im, c_mode_dtA_02_im);
        current_cell.store_vars(mode_dtA_03_re, c_mode_dtA_03_re);
        current_cell.store_vars(mode_dtA_03_im, c_mode_dtA_03_im);
        current_cell.store_vars(mode_dtA_04_re, c_mode_dtA_04_re);
        current_cell.store_vars(mode_dtA_04_im, c_mode_dtA_04_im);
        current_cell.store_vars(mode_dtA_05_re, c_mode_dtA_05_re);
        current_cell.store_vars(mode_dtA_05_im, c_mode_dtA_05_im);
        current_cell.store_vars(mode_dtA_06_re, c_mode_dtA_06_re);
        current_cell.store_vars(mode_dtA_06_im, c_mode_dtA_06_im);
        current_cell.store_vars(mode_dtA_07_re, c_mode_dtA_07_re);
        current_cell.store_vars(mode_dtA_07_im, c_mode_dtA_07_im);
        current_cell.store_vars(mode_dtA_08_re, c_mode_dtA_08_re);
        current_cell.store_vars(mode_dtA_08_im, c_mode_dtA_08_im);
        current_cell.store_vars(mode_dtA_09_re, c_mode_dtA_09_re);
        current_cell.store_vars(mode_dtA_09_im, c_mode_dtA_09_im);
        current_cell.store_vars(mode_dtA_10_re, c_mode_dtA_10_re);
        current_cell.store_vars(mode_dtA_10_im, c_mode_dtA_10_im);
        current_cell.store_vars(mode_dtA_11_re, c_mode_dtA_11_re);
        current_cell.store_vars(mode_dtA_11_im, c_mode_dtA_11_im);
        current_cell.store_vars(mode_dtA_12_re, c_mode_dtA_12_re);
        current_cell.store_vars(mode_dtA_12_im, c_mode_dtA_12_im);
        current_cell.store_vars(mode_dtA_13_re, c_mode_dtA_13_re);
        current_cell.store_vars(mode_dtA_13_im, c_mode_dtA_13_im);
        current_cell.store_vars(mode_dtA_14_re, c_mode_dtA_14_re);
        current_cell.store_vars(mode_dtA_14_im, c_mode_dtA_14_im);
        current_cell.store_vars(mode_dtA_15_re, c_mode_dtA_15_re);
        current_cell.store_vars(mode_dtA_15_im, c_mode_dtA_15_im);
        current_cell.store_vars(mode_dtA_16_re, c_mode_dtA_16_re);
        current_cell.store_vars(mode_dtA_16_im, c_mode_dtA_16_im);
        current_cell.store_vars(mode_dtA_17_re, c_mode_dtA_17_re);
        current_cell.store_vars(mode_dtA_17_im, c_mode_dtA_17_im);
        current_cell.store_vars(mode_dtA_18_re, c_mode_dtA_18_re);
        current_cell.store_vars(mode_dtA_18_im, c_mode_dtA_18_im);
        current_cell.store_vars(mode_dtA_19_re, c_mode_dtA_19_re);
        current_cell.store_vars(mode_dtA_19_im, c_mode_dtA_19_im);
        current_cell.store_vars(mode_dtA_20_re, c_mode_dtA_20_re);
        current_cell.store_vars(mode_dtA_20_im, c_mode_dtA_20_im);
    }
 };
  
 #endif /* MODEDECOMPOSITION_HPP */
 