//This file calculates the 3+1 decomposed components of the Stress Energy Tensor T_ab
//For a scalar field, with potential function as specified below
//It assume input in terms of the struct CCZ4Matter::vars_t defined in CCZ4Matter.hpp
#ifndef CCZ4EMTENSORSF_HPP_
#define CCZ4EMTENSORSF_HPP_

#include "TensorAlgebra.hpp"

template <class data_t>
struct emtensor_t
{
     tensor<2, data_t> Sij; // S_ij = T_ij
     tensor<1,data_t>  Si; // S_i = T_ia_n^a
     data_t            S; // S = S^i_i
     data_t            rho; // rho = T_ab n^a n^b
     data_t            dVdphi; // Gradient of V(phi)
};

class CCZ4EMTensorSF
{
    public:

    template <class data_t, template <typename> class vars_t>
    static emtensor_t<data_t>
    compute_emtensor_SF(
        const vars_t<data_t> &vars,
        const vars_t< tensor<1,data_t> >& d1,
        const tensor<2, data_t>& h_UU,
        const tensor<3, data_t>& chris,
        const data_t& advecphi
    )
    {
        emtensor_t<data_t> out;

    		// Calculate the stress energy tensor elements

        // Find the potential and its gradient in terms of phi
    		// do we want to specify in setup or always adjust here?
    		data_t Vofphi = vars.phi*vars.phi;		//WOULD LIKE TO BE ABLE TO USE SINE, COSINE, EXP ETC HERE 
    		out.dVdphi = 2.0*vars.phi;						//AND HERE

    		//components of stress energy tensor
    		data_t Vt = -vars.PiM * vars.PiM + 2.0*Vofphi;
    		FOR2(i,j)
    		{
        		Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    		}

    		data_t dphidt2 = (vars.lapse * vars.PiM + advecphi)*(vars.lapse * vars.PiM + advecphi);

    		FOR2(i,j)
				{ 
    				out.Sij[i][j] = -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
				}
    
    		out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij , h_UU);

    		tensor<1, data_t> T_i; // The T(i,3) components of the 4d stress energy tensor
    		FOR1(i)
    		{
     				T_i[i] = (d1.phi[i] * (vars.PiM*vars.lapse + advecphi));

      		FOR1(j)
      		{
   						T_i[i] += -0.5*Vt*vars.h[i][j]*vars.shift[j]/vars.chi;
      		}
    		}

    		FOR1(i)
    		{
     				out.Si[i] = - T_i[i]/vars.lapse;

      			FOR1(j)
      			{
   							out.Si[i] += vars.shift[j]/vars.lapse * out.Sij[i][j];
      			}
    		}

    		out.rho = dphidt2/vars.lapse/vars.lapse + 0.5*Vt; // = T_ab * n^a n^b
    		FOR2(i,j)
    		{
    			out.rho += (-0.5*Vt*vars.h[i][j]/vars.chi + out.Sij[i][j])*vars.shift[i]*vars.shift[j]/vars.lapse/vars.lapse;
    		}

        return out;
    }

};

#endif /* CCZ4EMTENSORSF_HPP_ */
