/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARSOLUTION_HPP_)
#error "This file should only be included through BosonStarSolution.hpp"
#endif

#ifndef BOSONSTARSOLUTION_IMPL_HPP_
#define BOSONSTARSOLUTION_IMPL_HPP_

// BosonStarSolution::BosonStarSolution() {}

inline BosonStarSolution::BosonStarSolution(BosonStar_params_t a_params_BosonStar,
                    Potential::params_t a_params_potential)
    : m_params_BosonStar(a_params_BosonStar), m_params_potential(a_params_potential)

{
}

void BosonStarSolution::main()
{
  	// finds the eigenvalue corresponding to intitial metric central values.
    // WW is big then it descends to the correct w+ and w- then uses interval
    // bisection to find best ww to machine precision

    // finds an eigenvalue with a lot of nodes, (20 + eigenstate) by defualt
  	WW = find_WW();
    // calculates the lower limit for eigenvalue to be used in interval bisection
  	lower_ww = ww_min(WW);
    // calculates the upper value for eigenvalue to be used in interval bisection
  	upper_ww = ww_max(WW,lower_ww);
    // does interval bisection
  	middle_ww = ww_IB(lower_ww,upper_ww);
  	ww = middle_ww;
    // finds the integer index at which the scalar field needs to be truncated
  	mid_int = find_midint();

    // force the scalar field to zero after the turning point and reintegrate the lapse and shift
    force_flat(mid_int);
    // (true) uses large radius adaptive stepsize to get asymptotics
    // (integrates vacuum metric out to huge radius ~ 10^10 to measure asymptotics).
    rk4_asymp(mid_int, true, ww);
    // now we have integrated to huge radius we know the asymptotics to high precision
    PSI_INF = psi[gridsize-1];
    OM_INF = omega[gridsize-1];


    // iterate process of rescaling radius 0 field values by the asymptotics
    // miraculously converges in about 2 iterations normally
    for (int q=0; q<5; q++)
    {
        PSC/=PSI_INF;
        OMC/=OM_INF;
        ww/=OM_INF*OM_INF;
        initialise();
        // the normal rk4 integration method (hidden inside all eigenvalue functions)
        rk4(ww);
        mid_int = find_midint();
        rk4_asymp(mid_int,true,ww);
        PSI_INF = psi[gridsize-1];
        OM_INF = omega[gridsize-1];
    }

    PSC/=PSI_INF;
    OMC/=OM_INF;
    ww/=OM_INF*OM_INF;
    initialise();
    rk4(ww);
    mid_int = find_midint();
    // (false) now do vacuum integration without adaptive stepsize for the
    // solution to be interpolated to grid
    rk4_asymp(mid_int,false,ww);

    // calculate the ADM mass and aspect mass at the edge of physical domain.
    // Aspect is more accurate but ADM should be similar.
    adm_mass = -psi[gridsize-1]*dpsi[gridsize-1]*radius_array[gridsize-1]*radius_array[gridsize-1];
    aspect_mass = 2.*radius_array[gridsize-1]*(sqrt(psi[gridsize-1])-1.);

    // optional print statement outputting star mass and eigenvalue
    std::cout << "Finished producing Bosonstar with ADM mass : " << adm_mass <<
    " aspect mass : "  << aspect_mass << " and eigenvalue w : " << sqrt(ww) <<
    std::endl;

}

// initialises the 5 field variables with their central values
void BosonStarSolution::initialise()
{	
  	p[0] = PC;
  	omega[0] = OMC;
  	psi[0] = PSC;
  	// field gradients must be zero at zero radius for physical solutions
  	dp[0] = 0.;
  	dpsi[0] = 0.;
    radius_array[0] = 0.;
}

// if the scalar field diverges it rounds down to a value wiht the same sign. This is to not affect axis crossings function (crossings) and let it accurately deal with inf/nan
void BosonStarSolution::fix()
{
    // borked loosely means broken
  	bool borked = false; // turns true if function gets over twice as large as central (r=0) value
  	double truncation;
  	for (int i = 0; i < gridsize; ++i)
  	{
    		if (not borked)
    		{
      			if (abs(p[i])> 2.0*p[0])
      			{
        				borked = true;
        				truncation = p[i];
      			}
    		}
    		if (borked)
    		{
    		    p[i] = truncation;
    		}
  	}
}

// sets scalar field and gradient to zero after the point decided by function find_midint
void BosonStarSolution::force_flat(const int iter_crit)
{
  	for (int i = iter_crit+1; i < gridsize; ++i)
  	{
    		p[i] = 0.;
    		dp[i] = 0.;
  	}
}


// finds the integer index at which the scalar field needs to be truncated
int BosonStarSolution::find_midint()
{
  	int crossings = 0, mid_int;

  	// follow the correct amount of crossings
  	for (int i = 0; i < gridsize-1; ++i)
  	{
    		if (crossings == EIGEN)
    		{
      			mid_int = i;
      			break;
    		}
    		if (p[i]*p[i+1] < 0.)
    		{
    			   crossings += 1;
    		}
  	}

  	// climb the hill if crossings != 0, this will exit loop immediately if there is no crossings
  	for (int i = mid_int+5; i < gridsize-1; ++i)
  	{
    		if (abs(p[i+1]) < abs(p[i]))
    		{
      			mid_int = i;
      			break;
    		}
  	}

  	for (int i = mid_int+5; i < gridsize-1; ++i)
  	{
    		if (abs(p[i+1]) > abs(p[i]))
    		{
      			mid_int = i;
      			//std::cout << "Truncation error: " << p[i]/PC << std::endl;
      			return mid_int;
    		}
  	}
  	return gridsize-1;
}


// finds an eigenvalue with a lot of nodes, (20 + eigenstate) by defualt
double BosonStarSolution::find_WW()
{
  	int eigenstate;
  	double WW_ = 1.;
  	while (true)
  	{
    		initialise();
    		rk4(WW_);
    		fix();
    		eigenstate = crossings();
    		if (eigenstate >= 20 + EIGEN){ return WW_;}
    		WW_*=2.;
  	}
}

// calculates the lower limit for eigenvalue to be used in interval bisection
double BosonStarSolution::ww_min(const double WW_)
{
  	int accuracy = 400, eigenstate;
  	double ww_, lower_ww_;
  	for (int i=0; i<accuracy; i++)
  	{
    		ww_ = WW_*(double)(accuracy - i)/(double)accuracy;
    		initialise();
    		rk4(ww_);
    		fix();
    		eigenstate = crossings();
    		if (eigenstate <= EIGEN)
    		{
      			lower_ww_ = ww_;
      			return lower_ww_;
    		}
  	}
  	return 0.;
}


// calculates the upper value for eigenvalue to be used in interval bisection
double BosonStarSolution::ww_max(const double WW_, const double lower_ww_)
{
  	int accuracy = 400, eigenstate;
  	double ww_, upper_ww_;
  	for (int i=0; i<accuracy; i++)
  	{
    		ww_ = lower_ww + WW_*((double)i)/((double)accuracy);
    		initialise();
    		rk4(ww_);
    		fix();
    		eigenstate = crossings();
    		if (eigenstate > EIGEN )
    		{
      			upper_ww_ = ww_;
      			return upper_ww_;
    		}
  	}
  	return 0.;
}
// calculates the upper value for eigenvalue to be used in interval bisection

// takes in an upper and lower eigenvalue and uses interval bisection to find the solution inbetween
double BosonStarSolution::ww_IB(double lower_ww_, double upper_ww_)
{
  	int iter = 0, itermax, eigenstate, decimal_places_of_omega = 25.;
  	double middle_ww_;


    //calculate number of bisections needed (simple pen and paper calculation)
  	itermax = (int)((log(upper_ww_)+decimal_places_of_omega*log(10.))/log(2.));
  	while (true)
  	{
    		iter++;
    		middle_ww_ = 0.5*(upper_ww_+lower_ww_);
    		initialise();
    		rk4(middle_ww_);
    		fix();
    		eigenstate = crossings();
    		if (eigenstate > EIGEN)
    		{
    			   upper_ww_ = middle_ww_;
    		}
    		else if (eigenstate <= EIGEN)
    		{
    			   lower_ww_ = middle_ww_;
    		}
    		if (iter > itermax){break;}
  	}
  	return middle_ww_;
}


// integrate the full ODE system from r=0 to dx*gridsize, this may (probably will) blow up at laarge raadius, but it is fixed by other functions aafterwards.
// it has smaller stepsize for the first (adaptive_buffer) steps.
void BosonStarSolution::rk4(const double ww_)
{
  	double k1, k2, k3, k4, q1, q2, q3, q4, x_=0., h = dx/2.;
  	const double DX = dx;
        double DX_ = DX;
  	double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4;
  	int index, jmax=0;
  	radius_array[0] = 0.;

  	for (int i = 1; i < gridsize; ++i)
    {
    		DX_ = DX;
    		jmax = 0;
    		if (i<adaptive_buffer)
    		{
    			   jmax = adaptive_stepsize_repetitions;
    		}
    		for (int j=0; j<=jmax; j++)
    		{
      			DX_ = DX/( (double)(1+jmax) );
      			h = DX_/2.;
            x_ = (i-1)*dx+j*DX_;

      			k1 = DX_*P_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
      			q1 = DX_*DP_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
      			o1 = DX_*OMEGA_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
      			s1 = DX_*PSI_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
      			r1 = DX_*DPSI_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);

      			k2 = DX_*P_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
      			q2 = DX_*DP_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
      	    o2 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
      			s2 = DX_*PSI_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
      			r2 = DX_*DPSI_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);

      			k3 = DX_*P_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
      			q3 = DX_*DP_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
            o3 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
      			s3 = DX_*PSI_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
      			r3 = DX_*DPSI_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);

      			k4 = DX_*P_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
      			q4 = DX_*DP_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
      	    o4 = DX_*OMEGA_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
      			s4 = DX_*PSI_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
      			r4 = DX_*DPSI_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);

      			index = i-1;
      			if(j==jmax)
      			{
      				index = i;
      			}

      			p[index] = p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
      			dp[index] = dp[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
      			psi[index] = psi[i-1] + (s1 + 2.*s2 + 2.*s3 + s4)/6.;
      			dpsi[index] = dpsi[i-1] + (r1 + 2.*r2 + 2.*r3 + r4)/6.;
      			omega[index] = omega[i-1] + (o1 + 2.*o2 + 2.*o3 + o4)/6.;
    		}
  	radius_array[i] = dx*i;
  	}

}

// takes an integrated ODE system and starts at point (iter) and re-integrates but enforcing scalara field to decaay or be in vacuum
// the integral is adaptive in that it aaccelerates ar later radius in order to find correct asymptotic behaviour. It will shout if the radius reached is below 10e10
// bool adaaptive is true if stepsize is supposed to be adaptive aat large radius and false for constant stepsize
void BosonStarSolution::rk4_asymp(const int iter, const bool adaptive, const double ww_)
{
    double k1, k2, k3, k4, q1, q2, q3, q4, x_=iter*dx, h, delta = (double)gridsize;
    const double DX = dx;
    double DX_= DX;
    double o1, o2, o3, o4, s1, s2, s3, s4, r1, r2, r3, r4;
    double N_ = gridsize-iter, L_ = pow(9.,9);
    int i_;

    double k_ = log(L_)/N_;

  	for (int i = iter+1; i < gridsize; ++i)
  	{
                i_ = double(i-iter);
    		if (adaptive)
    		{
      			if (x_<8e8)
      			{
      				    DX_ = (exp(k_)-1.)*exp(k_*i_);
      			}
      			else
      			{
      				    DX_ = DX;
      			}
    		}
    		h = DX_/2.;

    		//k1 = dx*small_P_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
    		//q1 = dx*DP_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
    		o1 = DX_*OMEGA_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
    		s1 = DX_*PSI_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);
    		r1 = DX_*DPSI_RHS(x_,p[i-1],dp[i-1],psi[i-1],dpsi[i-1],omega[i-1],ww_);

    		//k2 = dx*small_P_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
    		//q2 = dx*DP_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
        o2 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
    		s2 = DX_*PSI_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);
    		r2 = DX_*DPSI_RHS(x_ + h,p[i-1] + k1/2.,dp[i-1] + q1/2.,psi[i-1] + s1/2.,dpsi[i-1] + r1/2.,omega[i-1] + o1/2.,ww_);

    		//k3 = dx*small_P_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
    		//q3 = dx*DP_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
        o3 = DX_*OMEGA_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
    		s3 = DX_*PSI_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);
    		r3 = DX_*DPSI_RHS(x_ + h,p[i-1] + k2/2.,dp[i-1] + q2/2.,psi[i-1] + s2/2.,dpsi[i-1] + r2/2.,omega[i-1] + o2/2.,ww_);

    		//k4 = dx*small_P_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
    		//q4 = dx*DP_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
        o4 = DX_*OMEGA_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
    		s4 = DX_*PSI_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);
    		r4 = DX_*DPSI_RHS(x_ + 2.*h,p[i-1] + k3,dp[i-1] + q3,psi[i-1] + s3,dpsi[i-1] + r3,omega[i-1] + o3,ww_);

    		p[i] = 0.;//p[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
    		dp[i] = 0.;//dp[i-1] + (q1 + 2.*q2 + 2.*q3 + q4)/6.;
    		psi[i] = psi[i-1] + (s1 + 2.*s2 + 2.*s3 + s4)/6.;
    		dpsi[i] = dpsi[i-1] + (r1 + 2.*r2 + 2.*r3 + r4)/6.;
    		omega[i] = omega[i-1] + (o1 + 2.*o2 + 2.*o3 + o4)/6.;
    		x_ += DX_;
                if (!adaptive){radius_array[i]=i*dx;}
  	}

  	if (adaptive and x_ < 8e7)
  	{
  	        std::cout << "\33[30;41m" << " Asymptotic Radius Too Small" << "\x1B[0m" << std::endl;
                std::cout << x_ << std::endl;
  	}
}



// these functions return the right hand side of the ode's

// small_P_RHS is ONLY valid for large redius when the scalar field is small.

// can use this instead of vacuum at large radius for exponentially decaying scalar field. doesnt make much difference
double BosonStarSolution::small_P_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	double RHS = -P*PSI*sqrt(DV(P)-ww_/(OM*OM));
  	return RHS;
}

double BosonStarSolution::P_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	double RHS = DP;
  	return RHS;
}

double BosonStarSolution::DP_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	double DOM = OMEGA_RHS(x,P,DP,PSI,DPSI,OM,ww_), V = MM*P*P;
  	double RHS = P*PSI*PSI*(DV(P) - ww_/(OM*OM)) - DP*(DOM/OM + DPSI/PSI);
        if (x>10e-9)
        {
                 RHS += -2.*DP/x;
  	}
        return RHS;
}

double BosonStarSolution::PSI_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	double RHS = DPSI;
  	return RHS;
}

double BosonStarSolution::DPSI_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	double RHS = DPSI*DPSI/(2.*PSI) - 2*M_PI*G*(V(P)*pow(PSI,3) + DP*DP*PSI + ww_*PSI*pow(P*PSI/OM,2));
        if (x>10e-9)
        {
                RHS += - 2.*DPSI/x;
        }
  	return RHS;
}

double BosonStarSolution::OMEGA_RHS(const double x, const double P, const double DP, const double PSI, const double DPSI, const double OM, const double ww_)
{
  	return (4.*M_PI*G*x*pow(PSI*OM,2)*(DP*DP - PSI*PSI*V(P) + ww_*pow(P*PSI/OM,2)) - OM*OM*DPSI*(x*DPSI+2.*PSI) )/(2.*PSI*OM*(x*DPSI + PSI));
}


// V is klein gordon potential and DV is its gradient. Depends on star type in BosonStarParams.hpp
double BosonStarSolution::V(const double P)
{
  	if (!solitonic)
  	{
  		  return MM*P*P + 0.5* lambda*P*P*P*P;
  	}
  	else
  	{
  		  return MM*P*P*pow((1.-2.*pow(P/sigma,2)),2);
  	}
}
double BosonStarSolution::DV(const double P)
{
  	if (!solitonic)
  	{
  		  return MM + lambda*P*P;
  	}
  	else
  	{
  		  return MM-8.*MM*pow(P/sigma,2)+12.*MM*pow(P/sigma,4);
  	}
}

// counts how many times the function crosses the axis
int BosonStarSolution::crossings()
{
  	int number=0;
  	for (int i = 0; i < gridsize-1; ++i)
  	{
    		if (p[i]*p[i+1]<0)
    		{
    			   number += 1;
    		}
  	}
  	return number;
}

// 4th order error (cubic interpolation) for field. shouts if asked to fetch a value outside the ode solution
double BosonStarSolution::get_p_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?p[1]:p[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = p[iter];
    f3 = p[iter+1];
    f4 = p[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

double BosonStarSolution::get_dp_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?dp[1]:dp[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = dp[iter];
    f3 = dp[iter+1];
    f4 = dp[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

double BosonStarSolution::get_lapse_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?omega[1]:omega[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = omega[iter];
    f3 = omega[iter+1];
    f4 = omega[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}

double BosonStarSolution::get_psi_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?psi[1]:psi[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = psi[iter];
    f3 = psi[iter+1];
    f4 = psi[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}




double BosonStarSolution::get_dpsi_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?dpsi[1]:dpsi[iter-1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = dpsi[iter];
    f3 = dpsi[iter+1];
    f4 = dpsi[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline, from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./48.)*(f1 *(-3.+2.*a+12.*a*a-8.*a*a*a) +(3.+2.*a)*(-(1.+2.*a)*(-9.*f3+f4+6.*f3*a-2*f4*a)+3.*f2*(3.-8.*a+4.*a*a)));
    return interpolated_value;
}




double BosonStarSolution::get_dlapse_interp(const double r) const
{
    int iter = (int) floor(r/dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a = (r/dx)-floor(r/dx)-0.5; //fraction from midpoint of two values, a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter==0)?omega[1]:omega[iter-1]);
    f2 = omega[iter];
    f3 = omega[iter+1];
    f4 = omega[iter+2];

    if (iter>gridsize-3){std::cout << "Requested Value outside BS initial data domain!" << std::endl;}

    // do the cubic spline (for gradient now), from mathematica script written by Robin (rc634@cam.ac.uk)
    interpolated_value = (1./(24.*dx))*( (f1-27.*f2+27.*f3-f4)  +  12.*a*(f1-f2-f3+f4)  -  12.*a*a*(f1-3.*f2+3.*f3-f4)  );
    return interpolated_value;
}


// returns the aspect mass, i.e. comparing the large radius metric to Schwarzschild
// and using a differential realtion to return M (NOT an ADM mass calc)
double BosonStarSolution::get_mass() const
{
    return aspect_mass;
}

// returns the eigenvalue
double BosonStarSolution::get_w() const
{
    return sqrt(ww);
}


void BosonStarSolution::set_initialcondition_params(const double max_r)
{
    gridsize = m_params_BosonStar.gridpoints;
    adaptive_buffer = 0.;//gridsize/10; // numer of gridpoints to intergate more carefully
    p.resize(gridsize); //scalar field modulus
    dp.resize(gridsize); //scalar field modulus gradient
    psi.resize(gridsize); //conformal factor
    dpsi.resize(gridsize); //conformal factor gradient
    omega.resize(gridsize); // lapse
    radius_array.resize(gridsize); //radius

    G = m_params_BosonStar.Newtons_constant;
    PC = m_params_BosonStar.central_amplitude_CSF;
    EIGEN = m_params_BosonStar.eigen;
    MM = m_params_potential.scalar_mass*m_params_potential.scalar_mass;
    lambda = m_params_potential.phi4_coeff;
    solitonic = m_params_potential.solitonic;
    sigma = m_params_potential.sigma_soliton;
    L = max_r*1.05; //just to make sure the function domain is slightly larger than the required cube
    dx = L/(gridsize-1);
}


void BosonStarSolution::shout() const
{
    std::cout << "Haliboombah!" << std::endl;
}


#endif /* BOSONSTARSOLUTION_IMPL_HPP_ */
