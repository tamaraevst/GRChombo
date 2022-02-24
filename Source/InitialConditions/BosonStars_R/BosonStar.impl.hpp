/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_IMPL_HPP_
#define BOSONSTAR_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "WeightFunction.hpp"
#include "DebuggingTools.hpp"
#include "Max.hpp"

inline BosonStar::BosonStar(BosonStar_params_t a_params_BosonStar, BosonStar_params_t a_params_BosonStar2,
                    Potential::params_t a_params_potential, double a_G_Newton,
                    double a_dx, bool a_identical, int a_verbosity)
    :m_dx(a_dx), m_G_Newton(a_G_Newton), m_params_BosonStar(a_params_BosonStar), m_params_BosonStar2(a_params_BosonStar2),
    m_params_potential(a_params_potential), m_identical(a_identical), m_verbosity(a_verbosity)
{
}

void BosonStar::compute_1d_solution(const double max_r)
{
    try
    {   
        //set initial parameters and then run the solver (didnt put it in the constructor)
        
        pout() << "Setting initial conditions for Star 1" << endl;
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,m_params_potential,max_r);
        pout() << "Running the solver for Star 1" << endl;
        m_1d_sol.main();
        pout() << "Completed for star 1" << endl;

        pout() << "Setting initial conditions for Star 2" << endl;
        m_1d_sol2.set_initialcondition_params(m_params_BosonStar2,m_params_potential,max_r);
        pout() << "Running the solver for Star 2" << endl;
        m_1d_sol2.main();
        pout() << "Completed for star 2" << endl;

    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BosonStar::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);

    // Import BS parameters and option of whether this is a BS binary or BS-BH binary
    double rapidity = m_params_BosonStar.BS_rapidity;
    double rapidity2 = m_params_BosonStar2.BS_rapidity;
    bool binary = m_params_BosonStar.BS_binary;
    bool BS_BH_binary = m_params_BosonStar.BS_BH_binary;
    double M = m_params_BosonStar.BlackHoleMass;
    double separation = m_params_BosonStar.BS_separation;
    double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    double q = m_params_BosonStar.mass_ratio;
    double alpha = m_params_BosonStar.alpha_stretch;
    bool do_stretch = m_params_BosonStar.do_stretch;
    int n_weight = m_params_BosonStar.n_power;
    int initial_data_choice = m_params_BosonStar.id_choice;
    double radius_width1 = m_params_BosonStar.radius_width1;
    double radius_width2 = m_params_BosonStar.radius_width2;
    int conformal_power = m_params_BosonStar.conformal_factor_power;

    // Define boosts and coordinate objects, suppose star 1 is on the left of the centre of mass 
    // and star 2 is on the right of centre of mass

    //e.g. taking the centre of mass to be the origin, then STAR2 -------- (origin) -------- STAR1
    
    // First star positioning
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double v_ = tanh(rapidity);
    double t = (coords.x - separation / (q + 1.)) * s_; //set /tilde{t} to zero
    double x = (coords.x - separation / (q + 1.)) * c_;
    double z = coords.z; //set /tilde{t} to zero
    double y = coords.y + impact_parameter / 2.;
    double r = sqrt(x * x + y * y + z * z);

    // First star physical variables
    double p_ = m_1d_sol.get_p_interp(r);
    double dp_ = m_1d_sol.get_dp_interp(r);
    double omega_ = m_1d_sol.get_lapse_interp(r);
    double omega_prime_ = m_1d_sol.get_dlapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double psi_prime_ = m_1d_sol.get_dpsi_interp(r);

    // Get scalar field modulus, conformal factor, lapse and their gradients
    double pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
    double lapse_1 = omega_ * psi_ / (sqrt(pc_os));
    double lapse_2 = 1.;
    double w_ = m_1d_sol.get_w();

    //Write in phase, shift, metric componnets of star 1 and initialise metric components of star 2
    double phase_ = w_ * t;
    double beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    double g_zz_1 = psi_ * psi_;
    double g_yy_1 = psi_ * psi_;
    double g_xx_1 = pc_os;
    double g_xx_2 = 0., g_yy_2 = 0., g_zz_2 = 0., g_xx, g_yy, g_zz;

    //Add on to evolution equations
    vars.phi_Re += p_ * cos(phase_);
    vars.phi_Im += p_ * sin(phase_);
    vars.Pi_Re += -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) - w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
    vars.Pi_Im += -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) + w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

    //Initialise extrinsic curvature and metric with upper indices
    double KLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K1 = 0., K2 = 0.;

    // Fill them in
    gammaUU_1[0][0] = 1. / g_xx_1;
    gammaUU_1[1][1] = 1. / g_yy_1;
    gammaUU_1[2][2] = 1. / g_zz_1;

    KLL_1[2][2] = -lapse_1 * s_ * x * psi_prime_ / (r * psi_);
    KLL_1[1][1] = KLL_1[2][2];
    KLL_1[0][1] = lapse_1 * c_ * s_ * (y / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[0][2] = lapse_1 * c_ * s_ * (z / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = KLL_1[0][2];
    KLL_1[2][1] = 0.;
    KLL_1[1][2] = 0.;
    KLL_1[0][0] = lapse_1 * (x / r) * s_ * c_ * c_ * (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ + v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i,j) K1 += gammaUU_1[i][j] * KLL_1[i][j];

    // Here we use Thomas Helfer's trick and find the corresponding fixed values to be substracted in the initial guess
    double helferLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double helferLL2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    // Note that for equal mass helferLL = helferLL2
    
    // Here is the conformal factor that will be differently calculated depending on the choice of the initial data
    double chi_;
    double chi_plain;

     // This is the effect of object 1 on object 2 and hence represents the value to be substracted in the initial data from the position of object 2 
    double t_p = (-separation) * s_; //set /tilde{t} to zero
    double x_p = (-separation) * c_;
    double z_p = 0.; //set /tilde{t} to zero
    double y_p = impact_parameter;
    double r_p = sqrt(x_p * x_p + y_p * y_p + z_p * z_p);
    double p_p = m_1d_sol.get_p_interp(r_p);
    double dp_p = m_1d_sol.get_dp_interp(r_p);
    double omega_p = m_1d_sol.get_lapse_interp(r_p);
    double omega_prime_p = m_1d_sol.get_dlapse_interp(r_p);
    double psi_p = m_1d_sol.get_psi_interp(r_p);
    double psi_prime_p = m_1d_sol.get_dpsi_interp(r_p);
    double pc_os_p = psi_p * psi_p * c_ * c_ - omega_p * omega_p * s_ * s_;
    
    //Initialise weight function arguments to some random values - good check if in the compute
    //of weight functions these values should never appear
    // double arg1 = 42.0;
    // double arg2 = 42.0;

    // double stretch_factor1 = 1.0;
    // double stretch_factor2 = 1.0;

    WeightFunction weight;
    
    //double check_y = max(fabs(coords.y) - 2*separation, 0);
    //double check_z = max(fabs(coords.z) - 2*separation, 0);

    // if (do_stretch)
    // {
    //     double stretch_factor1 = weight.stretching_factor((coords.x - separation / (q+1)) * cosh(rapidity), coords.y, alpha);
    // }
    //Argument of weight function to be applied to star 1
	// arg1 = (stretch_factor1/separation) * (sqrt(pow((coords.x - separation / (q+1)) * cosh(rapidity), 2) + pow(coords.y,2) + pow(coords.z, 2)));

    if (binary)
    {
        helferLL[1][1] = psi_p * psi_p;
        helferLL[2][2] = psi_p * psi_p;
        helferLL[0][0] = pc_os_p;
        double chi_inf = pow((2. - helferLL[0][0]) * (2. - helferLL[1][1]) *
        (2. - helferLL[2][2]), -1./3.), h00_inf = (2. - helferLL[0][0]) * chi_inf,
        h11_inf = (2. - helferLL[1][1]) * chi_inf, h22_inf = (2. - helferLL[2][2]) * chi_inf;
        /*if (r<3){
        std::cout << "h00 = " << h00_inf << ", h11 = " << h11_inf
                          << ", h22 = " << h22_inf << ", chi inf = " <<
                          chi_inf << std::endl;}*/
    }

    if (binary)
    {   
        //Second star positioning
        c_ = cosh(-rapidity2);
        s_ = sinh(-rapidity2);
        v_ = tanh(-rapidity2);
        t = (coords.x + q * separation / (q + 1.)) * s_; //set /tilde{t} to zero
        x = (coords.x + q * separation / (q + 1.)) * c_;
        z = coords.z;
        y = coords.y - impact_parameter / 2.;
        r = sqrt(x * x + y * y + z * z);

        //Second star physical variables
        p_ = m_1d_sol2.get_p_interp(r);
        dp_ = m_1d_sol2.get_dp_interp(r);
        omega_ = m_1d_sol2.get_lapse_interp(r);
        omega_prime_ = m_1d_sol2.get_dlapse_interp(r);
        psi_ = m_1d_sol2.get_psi_interp(r);
        psi_prime_ = m_1d_sol2.get_dpsi_interp(r);
        double r_tilde;

        if (BS_BH_binary)
        {
            r_tilde = sqrt(r  *r + 10e-10);
            omega_ = (2. - M / r_tilde) / (2. + M / r_tilde);
            omega_prime_ = 4. * M / pow(2. * r_tilde + M, 2);
            psi_ = pow(1. + M/ (2. * r_tilde), 2);
            psi_prime_ = -(M / (r_tilde * r_tilde)) * (1. + M / (2. * r_tilde));
        }

        pc_os = psi_ * psi_ * c_ *c_ - omega_ * omega_ * s_ * s_;
        lapse_2 = omega_ * psi_ / (sqrt(pc_os));
        w_ = m_1d_sol2.get_w();
        phase_ = m_params_BosonStar.phase * M_PI + w_ * t;
        beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
        vars.shift[0] += beta_x;
        g_zz_2 = psi_ * psi_;
        g_yy_2 = psi_ * psi_;
        g_xx_2 = pc_os;
        gammaUU_2[0][0] = 1. / g_xx_2;
        gammaUU_2[1][1] = 1. / g_yy_2;
        gammaUU_2[2][2] = 1. / g_zz_2;

        if (BS_BH_binary)
        {
            //do not need to modify scalar field or momentum if we have a black hole
        }
        else
        {
            vars.phi_Re += p_ * cos(phase_);
            vars.phi_Im += p_ * sin(phase_);
            vars.Pi_Re += -(1. / lapse_2) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) - w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
            vars.Pi_Im += -(1. / lapse_2) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) + w_ * (c_ - beta_x * s_) * p_ * cos(phase_));
        }

        KLL_2[2][2] = -lapse_2 * s_ * x * psi_prime_ / (r * psi_);
        KLL_2[1][1] = KLL_2[2][2];
        KLL_2[0][1] = lapse_2 * c_ * s_ * (y / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
        KLL_2[0][2] = lapse_2 * c_ * s_ * (z / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
        KLL_2[1][0] = KLL_2[0][1];
        KLL_2[2][0] = KLL_2[0][2];
        KLL_2[2][1] = 0.;
        KLL_2[1][2] = 0.;
        KLL_2[0][0] = lapse_2 * (x / r) * s_ * c_ * c_ * (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ + v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
        FOR2(i,j) K2 += gammaUU_2[i][j] * KLL_2[i][j];

        // Again, finding the values to be substracted from position of star 1, that is below we find the effect of star 2 on star 1
        double x_p2 = (separation) * c_;
        double z_p2 = 0.; //set /tilde{t} to zero
        double y_p2 = -impact_parameter;
        double r_p2 = sqrt(x_p2 * x_p2 + y_p2 * y_p2 + z_p2 * z_p2);

        // Get physical variables needed for the metric
        double p_p2 = m_1d_sol2.get_p_interp(r_p2);
        double dp_p2 = m_1d_sol2.get_dp_interp(r_p2);
        double omega_p2 = m_1d_sol2.get_lapse_interp(r_p2);
        double omega_prime_p2 = m_1d_sol2.get_dlapse_interp(r_p2);
        double psi_p2 = m_1d_sol2.get_psi_interp(r_p2);
        double psi_prime_p2 = m_1d_sol2.get_dpsi_interp(r_p2);
        double pc_os_p2 = psi_p2 * psi_p2 * c_ * c_ - omega_p2 * omega_p2 * s_ * s_;

        if (m_identical == 1)
        {
            helferLL2[1][1] = psi_p * psi_p;
            helferLL2[2][2] = psi_p * psi_p;
            helferLL2[0][0] = pc_os_p;
        }
        else
        {
            helferLL2[1][1] = psi_p2 * psi_p2;
            helferLL2[2][2] = psi_p2 * psi_p2;
            helferLL2[0][0] = pc_os_p2;
        }
        
        // if (do_stretch)
        // {
        //     double stretch_factor2 = weight.stretching_factor2((coords.x + q * separation / (q + 1)) * cosh(-rapidity2), coords.y, alpha);
        // }
  
        //Argument of weight function to be applied to star 2
        // arg2 = (stretch_factor2/separation) * (sqrt(pow((coords.x + q * separation / (q + 1)) * cosh(-rapidity2), 2) + pow(coords.y,2) + pow(coords.z, 2)));
    }
    
    // Initial 3-metric 
    // g_xx = g_xx_1 + g_xx_2 - 1.0 - (weight1 * (helferLL[0][0] - 1.0) + weight2 * (helferLL2[0][0] - 1.0));
    // g_yy = g_yy_1 + g_yy_2 - 1.0 - (weight1 * (helferLL[1][1] - 1.0) + weight2 * (helferLL2[1][1] - 1.0));
    // g_zz = g_zz_1 + g_zz_2 - 1.0 - (weight1  * (helferLL[2][2] - 1.0) + weight2 * (helferLL2[2][2] - 1.0));

    // if (initial_data_choice !=3 )
    // {
    //     if (initial_data_choice == 0)
    //     {
    //         g_xx = g_xx_1 + g_xx_2 - 1.0;
    //         g_yy = g_yy_1 + g_yy_2 - 1.0;
    //         g_zz = g_zz_1 + g_zz_2 - 1.0;
    //     }

    //     if (initial_data_choice == 1)
    //     {
    //         g_xx = g_xx_1 + g_xx_2 - helferLL[0][0];
    //         g_yy = g_yy_1 + g_yy_2 - helferLL[1][1];
    //         g_zz = g_zz_1 + g_zz_2 - helferLL[2][2];
    //     }

    //     if (initial_data_choice == 5)
    //     {
    //         double weight1, weight2;

    //         // Use weight function for initial data. In case of BS-BH binary helferLL/helferLL2 varibales are zero so it doesn't make a difference there 
    
    //         weight1 = weight.compute_weight(arg1, n_weight); // bump at object 1
    //         weight2 = weight.compute_weight(arg2, n_weight); //bump at object 2

    //         //Just some sanity checks
    //         if (weight1 > 1.0)
    //         {DEBUG_OUT(weight1);}

    //         if (weight2 > 1.0)
    //         {DEBUG_OUT(weight2);}

    //         g_xx = g_xx_1 + g_xx_2 - 1.0 - (weight1 * (helferLL2[0][0] - 1.0) + weight2 * (helferLL[0][0] - 1.0));
    //         g_yy = g_yy_1 + g_yy_2 - 1.0 - (weight1 * (helferLL2[1][1] - 1.0) + weight2 * (helferLL[1][1] - 1.0));
    //         g_zz = g_zz_1 + g_zz_2 - 1.0 - (weight1  * (helferLL2[2][2] - 1.0) + weight2 * (helferLL[2][2] - 1.0));
    //     }

    //     if (initial_data_choice == 2)
    //     {
    //         double weight1, weight2;

    //         // Use weight function for initial data. In case of BS-BH binary helferLL/helferLL2 varibales are zero so it doesn't make a difference there 
    
    //         weight1 = weight.compute_weight(arg1, n_weight); // bump at object 1
    //         weight2 = weight.compute_weight(arg2, n_weight); //bump at object 2

    //         //Just some sanity checks
    //         if (weight1 > 1.0)
    //         {DEBUG_OUT(weight1);}

    //         if (weight2 > 1.0)
    //         {DEBUG_OUT(weight2);}

    //         double Htensor[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    //         double htensor[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

    //         Htensor[0][0] = (1.0/2.0) * (helferLL[0][0] + helferLL2[0][0] - 2.0);
    //         Htensor[1][1] = (1.0/2.0) * (helferLL[1][1] + helferLL2[1][1] - 2.0);
    //         Htensor[2][2] = (1.0/2.0) * (helferLL[2][2] + helferLL2[2][2] - 2.0);

    //         htensor[0][0] = (1.0/2.0) * (helferLL2[0][0] - helferLL[0][0]);
    //         htensor[1][1] = (1.0/2.0) * (helferLL2[1][1] - helferLL[1][1]);
    //         htensor[2][2] = (1.0/2.0) * (helferLL2[2][2] - helferLL[2][2]);

    //         g_xx = g_xx_1 + g_xx_2 - 1.0 - Htensor[0][0] - htensor[0][0] * (weight1 - weight2);
    //         g_yy = g_yy_1 + g_yy_2 - 1.0 - Htensor[1][1] - htensor[1][1] * (weight1 - weight2);
    //         g_zz = g_zz_1 + g_zz_2 - 1.0 - Htensor[2][2] - htensor[2][2] * (weight1 - weight2);
    //     }

    //     // Now, compute upper and lower components
    //     gammaLL[0][0] = g_xx;
    //     gammaLL[1][1] = g_yy;
    //     gammaLL[2][2] = g_zz;
    //     gammaUU[0][0] = 1./g_xx;
    //     gammaUU[1][1] = 1./g_yy;
    //     gammaUU[2][2] = 1./g_zz;

    //     // Define initial conformal factor
    //     chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

    //     chi_plain = chi_;
    // }
    
   
    //If one uses fixing conformal trick, we need to have the vales of the metric of star 1 at its centre
    double r_11 = 0.;
    double p_11 = m_1d_sol.get_p_interp(r_11);
    double dp_11 = m_1d_sol.get_dp_interp(r_11);
    double omega_11 = m_1d_sol.get_lapse_interp(r_11);
    double omega_prime_11 = m_1d_sol.get_dlapse_interp(r_11);
    double psi_11 = m_1d_sol.get_psi_interp(r_11);
    double psi_prime_11 = m_1d_sol.get_dpsi_interp(r_11);
    double pc_os_11 = psi_11 * psi_11 * cosh(rapidity) * cosh(rapidity) - omega_11 * omega_11 * sinh(rapidity) * sinh(rapidity);

    //If one uses fixing conformal trick, we need to have the vales of the metric of star 2 at its centre
    double r_22 = 0.;
    double p_22 = m_1d_sol2.get_p_interp(r_22);
    double dp_22 = m_1d_sol2.get_dp_interp(r_22);
    double omega_22 = m_1d_sol2.get_lapse_interp(r_22);
    double omega_prime_22 = m_1d_sol2.get_dlapse_interp(r_22);
    double psi_22 = m_1d_sol2.get_psi_interp(r_22);
    double psi_prime_22 = m_1d_sol2.get_dpsi_interp(r_22);
    double pc_os_22 = psi_22 * psi_22 * cosh(-rapidity2) * cosh(-rapidity2) - omega_22 * omega_22 * sinh(-rapidity2) * sinh(-rapidity2);

    //These are to be filled in with plain supporposed metric components evaluated at x_A and x_B respectively 
    double superpose_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double superpose_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
        
    //metric components of \gamma_A(x_A)
    double g_zz_11 = psi_11 * psi_11;
    double g_yy_11 = psi_11 * psi_11;
    double g_xx_11 = pc_os_11;

    //metric components of \gamma_B(x_B)
    double g_zz_22 = psi_22 * psi_22;
    double g_yy_22 = psi_22 * psi_22;
    double g_xx_22 = pc_os_22;

    // This  is \gamma_{ij}(x_A) = \gamma_A(x_A) + \gamma_B(x_A) - 1
    superpose_1[0][0] = g_xx_11 + helferLL2[0][0] - 1.;
    superpose_1[1][1] = g_yy_11 + helferLL2[1][1] - 1.;
    superpose_1[2][2] = g_zz_11 + helferLL2[2][2] - 1.;

    // This  is \gamma_{ij}(x_B) = \gamma_B(x_B) + \gamma_A(x_B) - 1
    superpose_2[0][0] = g_xx_22 + helferLL[0][0] - 1.;
    superpose_2[1][1] = g_yy_22 + helferLL[1][1] - 1.;
    superpose_2[2][2] = g_zz_22 + helferLL[2][2] - 1.;

    double n_power = conformal_power / 12.0;

    //This is \chi(x_A)
    double chi_1 = pow(superpose_1[0][0] * superpose_1[1][1] * superpose_1[2][2], n_power);
    //This is \chi(x_B)
    double chi_2 = pow(superpose_2[0][0] * superpose_2[1][1] * superpose_2[2][2], n_power);

    //This is \chi^A(x_A)
    double chi1_1 = pow(g_xx_11 * g_yy_11 * g_zz_11, n_power);
    //This is \chi^B(x_B)
    double chi2_2 = pow(g_xx_22 * g_yy_22 * g_zz_22, n_power);
        
    //This is \delta_A
    double delta_1 = chi1_1 - chi_1;
    //This is \delta_B
    double delta_2 = chi2_2 - chi_2;

    //Plain superposition 
    g_xx = g_xx_1 + g_xx_2 - 1.0;
    g_yy = g_yy_1 + g_yy_2 - 1.0;
    g_zz = g_zz_1 + g_zz_2 - 1.0;

    chi_plain = pow(g_xx * g_yy * g_zz, n_power);

    //Find all the profile functions needed
    double profile1 = weight.profile_chi((coords.x - separation / (q+1)) * cosh(rapidity), coords.y, coords.z, radius_width1);
    double profile2 = weight.profile_chi((coords.x + q * separation / (q+1)) * cosh(-rapidity2), coords.y, coords.z, radius_width2);
        
    double profile_11 = weight.profile_chi(0., 0., 0., radius_width1);
    double argument_xB_xA = (separation / (q+1)) * (-q * cosh(-rapidity2) - cosh(rapidity));
    double profile_12 = weight.profile_chi(argument_xB_xA, 0., 0., radius_width1);
        
    double argument_xA_xB = (separation / (q + 1)) * (cosh(rapidity) + q * cosh(-rapidity2));
    double profile_21 = weight.profile_chi(argument_xA_xB, 0., 0., radius_width2);
    double profile_22 = weight.profile_chi(0., 0., 0., radius_width2);

    double value1 = (-profile_21 * delta_2 + profile_22 * delta_1)/(profile_11 * profile_22 - profile_12 * profile_21);
    double value2 = (profile_11 * delta_2 - profile_12 * delta_1)/(profile_11 * profile_22 - profile_12 * profile_21);

    chi_ = chi_plain + profile1 * value1 + profile2 * value2;

    // Now, compute upper and lower components
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    // vars.chi = chi_;

    // g_xx = g_xx_1 + g_xx_2 - 1.0;
    // g_yy = g_yy_1 + g_yy_2 - 1.0;
    // g_zz = g_zz_1 + g_zz_2 - 1.0;

    //Now, compute upper and lower components
    // gammaLL[0][0] = g_xx;
    // gammaLL[1][1] = g_yy;
    // gammaLL[2][2] = g_zz;
    // gammaUU[0][0] = 1. / g_xx;
    // gammaUU[1][1] = 1. / g_yy;
    // gammaUU[2][2] = 1. / g_zz;

    // Define initial conformal factor
    // chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

    // chi_plain = chi_;

    double factor = pow((chi_plain / chi_), 4. / n_power);

    vars.chi = pow(chi_, - 4.0 / conformal_power);

    // Define initial lapse
    if (BS_BH_binary){vars.lapse += sqrt(vars.chi);}
    else if (binary){vars.lapse += sqrt(lapse_1 * lapse_1 + lapse_2 * lapse_2 - 1.);}
    else{vars.lapse += lapse_1;}

    // Define initial trace of K and A_ij
    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = pow(chi_plain, - 4.0 / conformal_power) * gammaLL[i][j];
    FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] + gammaUU_2[l][k] * KLL_2[k][j]);
    FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
    FOR2(i,j) vars.A[i][j] = pow(chi_plain, - 4.0 / conformal_power)  * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

    current_cell.store_vars(factor, c_factor);
    current_cell.store_vars(vars);
}

#endif /* BOSONSTAR_IMPL_HPP_ */