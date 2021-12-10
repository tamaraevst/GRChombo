/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBS_HPP_)
#error "This file should only be included through BinaryBS.hpp"
#endif

#ifndef BINARYBS_IMPL_HPP_
#define BINARYBS_IMPL_HPP_

#include "DebuggingTools.hpp"

inline BinaryBS::BinaryBS(BosonStar_params_t a_bosonstar_params, BosonStar_params_t a_bosonstar2_params,
                    Potential::params_t a_params_potential, double a_G_Newton,
                    double a_dx, int a_verbosity)
    :m_dx(a_dx), m_G_Newton(a_G_Newton), m_bosonstar(a_bosonstar_params, a_params_potential),m_bosonstar2(a_bosonstar2_params, a_params_potential),
    m_params_potential(a_params_potential), m_verbosity(a_verbosity)
{
}

void BinaryBS::compute_1d_solution(const double max_r)
{
    try
    {
        // Set initial parameters and then run the solver for both of the BSs (didnt put it in the constructor)
        m_bosonstar.set_initialcondition_params(max_r);
        m_bosonstar.main();

        m_bosonstar2.set_initialcondition_params(max_r);
        m_bosonstar2.main();
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BinaryBS::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)

    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_bosonstar.m_params_BosonStar.star_centre);

    // Import BS parameters andd option of whether this is a BS binary or BS-BH binary
    double rapidity = m_bosonstar.m_params_BosonStar.BS_rapidity;
    std::cout << "Rapidity 1 is" << rapidity << std::endl;
    double rapidity2 = m_bosonstar2.m_params_BosonStar.BS_rapidity;
    std::cout << "Rapidity 2 is" << rapidity2 << std::endl;
    double mu = m_bosonstar.m_params_BosonStar.mass_ratio;
    double M = m_bosonstar.m_params_BosonStar.BlackHoleMass;
    double separation = m_bosonstar.m_params_BosonStar.BS_separation;
    double impact_parameter = m_bosonstar.m_params_BosonStar.BS_impact_parameter;
    bool BS_binary = m_bosonstar.m_params_BosonStar.BS_binary;
    bool BS_BH_binary = m_bosonstar.m_params_BosonStar.BS_BH_binary;

    // Define boosts and coordinate objects

    // First star positioning 
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double v_ = tanh(rapidity);
    double t1 = (coords.x-separation)*s_; //set /tilde{t} to zero
    double x1 = (coords.x-separation)*c_;
    double z1 = coords.z; //set /tilde{t} to zero
    double y1 = coords.y+impact_parameter/2.;
    double r1 = sqrt(x1*x1+y1*y1+z1*z1);

    // Second star positioning 
    double c_2 = cosh(-rapidity2);
    double s_2 = sinh(-rapidity2);
    double v_2 = tanh(-rapidity2);
    double t2 = (coords.x+separation/2.)*s_2; //set /tilde{t} to zero
    double x2 = (coords.x+separation/2.)*c_2;
    double z2 = coords.z;
    double y2 = coords.y-impact_parameter/2.;
    double r2 = sqrt(x2*x2+y2*y2+z2*z2);

    // First star physical variables

    // Get scalar field modulus, conformal factor, lapse and their gradients
    double p_ = m_bosonstar.get_p_interp(r1);
    double dp_ = m_bosonstar.get_dp_interp(r1);
    double omega_ = m_bosonstar.get_lapse_interp(r1);
    double omega_prime_ = m_bosonstar.get_dlapse_interp(r1);
    double psi_ = m_bosonstar.get_psi_interp(r1);
    double psi_prime_ = m_bosonstar.get_dpsi_interp(r1);

    // Get lapse and shift
    double pc_os = psi_*psi_*c_*c_ - omega_*omega_*s_*s_; //this is denominator for the lapse in boosted frame and also 11 component of the 3-metric
    double lapse_1 = omega_*psi_/(sqrt(pc_os));
    double lapse_2 = 1.;
    double beta_x = s_*c_*(psi_*psi_-omega_*omega_)/(pc_os);
    vars.shift[0] += beta_x;

    // Get frequency and phase for star 1
    double w_ = m_bosonstar.get_w();
    double phase_ = w_*t1;

    // Get 3-metric components for the fist star in boosted frame and put them into the matrix
    double g_zz_1 = psi_*psi_; 
    double g_yy_1 = psi_*psi_;
    double g_xx_1 = pc_os;

    double gammaUU_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    gammaUU_1[0][0] = 1./g_xx_1;
    gammaUU_1[1][1] = 1./g_yy_1;
    gammaUU_1[2][2] = 1./g_zz_1;

    // Define the total metric and second object's metric 
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //this is total metric with lower indices
    double gammaUU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //this is total metric with upper indices
    double gammaUU_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //this one will be filled in later

    // Componets of the second metric
    double g_xx_2=0., g_yy_2=0., g_zz_2=0.; 

    // Components of the total metric
    double g_xx, g_yy, g_zz;

    // Evolution equations for Re and Im parts of scalar field and conjugate momentum in boosted frame
    vars.phi_Re += p_*cos(phase_);
    vars.phi_Im += p_*sin(phase_);
    vars.Pi_Re += -(1./lapse_1)*( (x1/r1)*(s_-beta_x*c_)*dp_*cos(phase_) - w_*(c_-beta_x*s_)*p_*sin(phase_) );
    vars.Pi_Im += -(1./lapse_1)*( (x1/r1)*(s_-beta_x*c_)*dp_*sin(phase_) + w_*(c_-beta_x*s_)*p_*cos(phase_) );

    // Define extrinsic curvature for both of the objects and the total one
    double KLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //for object 1
    double KLL_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //for objeect 2
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}; //total
    KLL_1[2][2] = -lapse_1*s_*x1*psi_prime_/(r1*psi_);
    KLL_1[1][1] = KLL_1[2][2];
    KLL_1[0][1] = lapse_1*c_*s_*(y1/r1)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL_1[0][2] = lapse_1*c_*s_*(z1/r1)*(psi_prime_/psi_ - omega_prime_/omega_ );
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = KLL_1[0][2];
    KLL_1[2][1] = 0.;
    KLL_1[1][2] = 0.;
    KLL_1[0][0] = lapse_1*(x1/r1)*s_*c_*c_*(psi_prime_/psi_ - 2.*omega_prime_/omega_ + v_*v_*omega_*omega_prime_*pow(psi_,-2));

    double K1=0., K2=0.; //these are componnet of K with upper indices
    FOR2(i,j) K1 += gammaUU_1[i][j]*KLL_1[i][j];

    // Here we use Thomas Helfer's trick and find the corresponding fixed values to be substracted in the initial guess
    double helferLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double helferLL2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

    // Note that for equal mass helferLL = helferLL2

    if (BS_binary) //if BS binary or BS/BH binary
    {
        // This is the effect of object 2 on object 1 and hence represents the value to be substracted in the initial data
        double t_p = -(separation)*s_; //set /tilde{t} to zero
        double x_p = -(separation)*c_;
        double z_p = 0.; //set /tilde{t} to zero
        double y_p = impact_parameter;
        double r_p = sqrt(x_p*x_p+y_p*y_p+z_p*z_p);
    
        // Run BS solver to be able to find lapse and conformal factor needed for the metric 
        double p_p = m_bosonstar.get_p_interp(r_p);
        double dp_p = m_bosonstar.get_dp_interp(r_p);
        double omega_p = m_bosonstar.get_lapse_interp(r_p);
        double omega_prime_p = m_bosonstar.get_dlapse_interp(r_p);
        double psi_p = m_bosonstar.get_psi_interp(r_p);
        double psi_prime_p = m_bosonstar.get_dpsi_interp(r_p);
        double pc_os_p = psi_p*psi_p*c_*c_ - omega_p*omega_p*s_*s_;

        //These are the values to be substracted in the initial data
        helferLL[1][1] = psi_p*psi_p;
        helferLL[2][2] = psi_p*psi_p;
        helferLL[0][0] = pc_os_p;

        // second star physical variables
        double p_2 = m_bosonstar2.get_p_interp(r2);
        double dp_2 = m_bosonstar2.get_dp_interp(r2);
        double omega_2 = m_bosonstar2.get_lapse_interp(r2);
        double omega_prime_2 = m_bosonstar2.get_dlapse_interp(r2);
        double psi_2 = m_bosonstar2.get_psi_interp(r2);
        double psi_prime_2 = m_bosonstar2.get_dpsi_interp(r2);
        // double r_tilde;

        // Again, finding the values to be substracted from this second star
        double t_p2 = (separation)*s_2; //set /tilde{t} to zero
        double x_p2 = (separation)*c_2;
        double z_p2 = 0.; //set /tilde{t} to zero
        double y_p2 = -impact_parameter;
        double r_p2 = sqrt(x_p2*x_p2+y_p2*y_p2+z_p2*z_p2);

        // Get physical variables needed for the metric
        double p_p2 = m_bosonstar2.get_p_interp(r_p2);
        double dp_p2 = m_bosonstar2.get_dp_interp(r_p2);
        double omega_p2 = m_bosonstar2.get_lapse_interp(r_p2);
        double omega_prime_p2 = m_bosonstar2.get_dlapse_interp(r_p2);
        double psi_p2 = m_bosonstar2.get_psi_interp(r_p2);
        double psi_prime_p2 = m_bosonstar2.get_dpsi_interp(r_p2);
        double pc_os_p2 = psi_p2*psi_p2*c_2*c_2 - omega_p2*omega_p2*s_2*s_2;

        // Values to be substracted 
        helferLL2[1][1] = psi_p2*psi_p2;
        helferLL2[2][2] = psi_p2*psi_p2;
        helferLL2[0][0] = pc_os_p2;

        if (BS_BH_binary)
        {
            double r_tilde = sqrt(r2*r2 + 10e-10);
            omega_2 = (2.-M/r_tilde)/(2.+M/r_tilde);
            omega_prime_2 = 4.*M/pow(2.*r_tilde + M,2);
            psi_2 = pow(1.+M/(2.*r_tilde),2);
            psi_prime_2 = -(M/(r_tilde*r_tilde))*(1.+M/(2.*r_tilde));
        }

        // Find lapse, shift for object 2
        double pc_os2 = psi_2*psi_2*c_2*c_2 - omega_2*omega_2*s_2*s_2;
        lapse_2 = omega_2*psi_2/(sqrt(pc_os2));
        double beta_x2 = s_2*c_2*(psi_2*psi_2-omega_2*omega_2)/(pc_os2);
        vars.shift[0] += beta_x2;
 
        // Find frequency and phase for object 2
        double w_2 = m_bosonstar2.get_w();
        double phase_2 = m_bosonstar2.m_params_BosonStar.phase*M_PI + w_2*t2;

        // Metric components for object 2
        g_zz_2 = psi_2*psi_2;
        g_yy_2 = psi_2*psi_2;
        g_xx_2 = pc_os2;
        gammaUU_2[0][0] = 1./g_xx_2;
        gammaUU_2[1][1] = 1./g_yy_2;
        gammaUU_2[2][2] = 1./g_zz_2;


        if (BS_BH_binary)
        {
            //do not need to modify scalar field or momentum if we have a black hole
        }
        else
        {
            // For BS binary need to add on the second star's variables to the evolution equations
            vars.phi_Re += p_2*cos(phase_2);
            vars.phi_Im += p_2*sin(phase_2);
            vars.Pi_Re += -(1./lapse_2)*( (x2/r2)*(s_2*(-beta_x2)*c_2)*dp_2*cos(phase_2) - w_2*(c_2*(-beta_x2)*s_2)*p_2*sin(phase_2) );
            vars.Pi_Im += -(1./lapse_2)*( (x2/r2)*(s_2*(-beta_x2)*c_2)*dp_2*sin(phase_2) + w_2*(c_2*(-beta_x2)*s_2)*p_2*cos(phase_2) );
        }

        // Find extrinsic curvature for object 2
        KLL_2[2][2] = -lapse_2*s_2*x2*psi_prime_2/(r2*psi_2);
        KLL_2[1][1] = KLL_2[2][2];
        KLL_2[0][1] = lapse_2*c_2*s_2*(y2/r2)*(psi_prime_2/psi_2 - omega_prime_2/omega_2 );
        KLL_2[0][2] = lapse_2*c_2*s_2*(z2/r2)*(psi_prime_2/psi_2 - omega_prime_2/omega_2 );
        KLL_2[1][0] = KLL_2[0][1];
        KLL_2[2][0] = KLL_2[0][2];
        KLL_2[2][1] = 0.;
        KLL_2[1][2] = 0.;
        KLL_2[0][0] = lapse_2*(x2/r2)*s_2*c_2*c_2*(psi_prime_2/psi_2 - 2.*omega_prime_2/omega_2 + v_2*v_2*omega_2*omega_prime_2*pow(psi_2,-2));
        FOR2(i,j) K2 += gammaUU_2[i][j]*KLL_2[i][j];

    }
    
    // Use weight function for initial data. In case of BS-BH binary helferLL/helferLL2 varibales are zero so it doesn't make a difference there 
    WeightFunction weight;

    double arg1 = (1/separation) * (sqrt(pow(coords.x-x1, 2)+pow(coords.y-y1,2)));
    double arg2 = (1/separation) * (sqrt(pow(coords.x-x2, 2)+pow(coords.y-y2,2)));

    double weight1 = weight.weightfunction(arg1); // bump at object 1
    double weight2 = weight.weightfunction(arg2); //bump at object 2

    // Initial 3-metric 
    g_xx = g_xx_1 + g_xx_2 - 1.0 - (weight1 * (helferLL[0][0] - 1.0) + weight2 * (helferLL2[0][0] - 1.0));
    g_yy = g_yy_1 + g_yy_2 - 1.0 - (weight1 * (helferLL[1][1] - 1.0) + weight2 * (helferLL2[1][1] - 1.0));
    g_zz = g_zz_1 + g_zz_2 - 1.0 - (weight1 * (helferLL[2][2] - 1.0) + weight2 * (helferLL2[2][2] - 1.0));

    //These are the asymptotics that are used for Helfer trick, in our case asymptotoc values of chi and diagonal h are exacly 1 by construction
    // double chi_inf = pow((2.-helferLL[0][0])*(2.-helferLL[1][1])*
    //     (2.-helferLL[2][2]),-1./3.), h00_inf = (2.-helferLL[0][0])*chi_inf,
    //     h11_inf = (2.-helferLL[1][1])*chi_inf, h22_inf = (2.-helferLL[2][2])*chi_inf;
        /*if (r<3){
        std::cout << "h00 = " << h00_inf << ", h11 = " << h11_inf
                          << ", h22 = " << h22_inf << ", chi inf = " <<
                          chi_inf << std::endl;}*/

    // Now, compute upper and lower components
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1./g_xx;
    gammaUU[1][1] = 1./g_yy;
    gammaUU[2][2] = 1./g_zz;

    // Define initial conformal factor 
    double chi_ = pow(g_xx*g_yy*g_zz,-1./3.);
    vars.chi = chi_;

    // Define initial lapse 
    if (BS_BH_binary){vars.lapse += sqrt(vars.chi);}
    else if (BS_binary){vars.lapse += sqrt(lapse_1*lapse_1 + lapse_2*lapse_2-1.);}
    else{vars.lapse += lapse_1;}

    // Define initial trace of K and A_ij
    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = vars.chi*gammaLL[i][j];
    FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l]*(gammaUU_1[l][k]*KLL_1[k][j] + gammaUU_2[l][k]*KLL_2[k][j]);
    FOR2(i,j) vars.K += KLL[i][j]*gammaUU[i][j];
    FOR2(i,j) vars.A[i][j] = chi_*(KLL[i][j]-one_third*vars.K*gammaLL[i][j]);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
    current_cell.store_vars(weight1, c_weight1);
    current_cell.store_vars(weight2, c_weight2);
}

#endif /* BINARYBS_IMPL_HPP_ */
