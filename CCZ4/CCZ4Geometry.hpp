//This file calculates CCZ4 geometric quantities (or a similar 3+1 split.
//It assume input in terms of the struct vars_t defined in CCZ4.hpp

#include "CCZ4.hpp"

template <class data_t>
tensor<2, data_t>
compute_inverse_metric(const vars_t<data_t> &vars)
{
    data_t deth = vars.h[0][0]*vars.h[1][1]*vars.h[2][2] + 2*vars.h[0][1]*vars.h[0][2]*vars.h[1][2] - vars.h[0][0]*vars.h[1][2]*vars.h[1][2] - vars.h[1][1]*vars.h[0][2]*vars.h[0][2] - vars.h[2][2]*vars.h[0][1]*vars.h[0][1];
    tensor<2, data_t> h_UU;
    {
        h_UU[0][0] = (vars.h[1][1]*vars.h[2][2] - vars.h[1][2]*vars.h[1][2]) / deth;
        h_UU[0][1] = (vars.h[0][2]*vars.h[1][2] - vars.h[0][1]*vars.h[2][2]) / deth;
        h_UU[0][2] = (vars.h[0][1]*vars.h[1][2] - vars.h[0][2]*vars.h[1][1]) / deth;
        h_UU[1][1] = (vars.h[0][0]*vars.h[2][2] - vars.h[0][2]*vars.h[0][2]) / deth;
        h_UU[1][2] = (vars.h[0][1]*vars.h[0][2] - vars.h[0][0]*vars.h[1][2]) / deth;
        h_UU[2][2] = (vars.h[0][0]*vars.h[1][1] - vars.h[0][1]*vars.h[0][1]) / deth;
        h_UU[1][0] = h_UU[0][1];
        h_UU[2][0] = h_UU[0][2];
        h_UU[2][1] = h_UU[1][2];
    }
    return h_UU;
}

struct chris_t
{
    tensor<3, data_t> LLL; //3 lower indices
    tensor<3, data_t> ULL;
    tensor<1, data_t> contracted;  //chris contracted

    template <class data_t>
    chris_t(const vars_t<data_t> &vars,const vars_t<data_t> (&d1)[CH_SPACEDIM],const tensor<2, data_t>& h_UU)
    {
        FOR3(i,j,k)
        {
            LLL[i][j][k] = 0.5*(d1[k].h[j][i] + d1[j].h[k][i] - d1[i].h[j][k]);
        }
        FOR3(i,j,k)
        {
            ULL[i][j][k] = 0;
            FOR1(l)
            {
                ULL[i][j][k] += h_UU[i][l]*LLL[l][j][k];
            }
        }

    // Technically we can write contracted[i] = h_UU[j][k]*chris[i][j][k],
    // but this is not numerically stable: h_UU[j][k]*d1[i].h[j][k] should be zero
    // but in practice can be > O(1).
        FOR1(i)
        {
            contracted[i] = 0;
            FOR3(j,k,l)
            {
                contracted[i] += h_UU[i][j]*h_UU[k][l]*d1[l].h[k][j];
            }
        }
    }

    chris_t(const vars_t<data_t> &vars,const vars_t<data_t> (&d1)[CH_SPACEDIM])
    {
      const tensor<2, data_t> h_UU = compute_inverse_metric(vars);
      chris_t(vars,d1,h_UU);
    }
};

////more readable in the main code
//template <class data_t>
//christ_t
//compute_christoffel(const vars_t<data_t> &vars,const vars_t<data_t> (&d1)[CH_SPACEDIM],const tensor<2, data_t>& h_UU)
//{
//   return chris_t(vars,d1,h_UU);
//}
//
//template <class data_t>
//christ_t
//compute_christoffel(const vars_t<data_t> &vars,const vars_t<data_t> (&d1)[CH_SPACEDIM])
//{
//   const tensor<2, data_t> h_UU = compute_inverse_metric(vars);
//   return chris_t(vars,d1,h_UU);
//}

// Trick: For CCZ4, we can add Z terms to ricci by changing Gamma to chrisvec
template <class data_t>
struct ricciZ_t
{
   tensor<2, data_t> LL; //Ricci with two indices down
   data_t            scalar; //Ricci scalar

   ricciZ_t(const vars_t<data_t> &vars, const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
         const vars_t<data_t> (&d1)[CH_SPACEDIM],
         const christ_t& chris,
         const tensor<2, data_t>& h_UU,
         const tensor<1, data_t>& Z_over_chi)
   {
      data_t boxtildechi = 0;
      FOR2(k,l)
      {
         boxtildechi += h_UU[k][l]*covdtilde2chi[k][l];
      }

      tensor<2, data_t> covdtilde2chi;
      {
         FOR2(k,l)
         {
            covdtilde2chi[k][l] = d2[k][l].chi;
            FOR1(m)
            {
               covdtilde2chi[k][l] -= chris[m][k][l]*d1[m].chi;
            }
         }
      }

      FOR2(i,j)
      {
         data_t ricci_tilde = 0;
         FOR1(k)
         {
            ricci_tilde += 0.5*(vars.h[k][i]*d1[j].Gamma[k] + vars.h[k][j]*d1[i].Gamma[k];
            //This way of writing it allows the user to pass Z/chi = {0};
            ricci_tilde += (vars.Gamma[k] - 2*Z_over_chi[k]) chrisvec[k]*(chris_LLL[i][j][k] + chris_LLL[j][i][k]));
            FOR1(l)
            {
               ricci_tilde -= 0.5*h_UU[k][l]*d2[k][l].h[i][j];
               FOR1(m)
               {
                  ricci_tilde += h_UU[l][m]*(chris[k][l][i]*chris_LLL[j][k][m] + chris[k][l][j]*chris_LLL[i][k][m] + chris[k][i][m]*chris_LLL[k][l][j]);
               }
            }
         }

         data_t ricci_chi = 0.5*((GR_SPACEDIM-2)*covdtilde2chi[i][j] + vars.h[i][j]*boxtildechi - ((GR_SPACEDIM-2)*d1[i].chi*d1[j].chi + GR_SPACEDIM*vars.h[i][j]*dchi_dot_dchi) / (2*chi_regularised));

         data_t z_terms = 0;
         FOR1(k)
         {
            z_terms += Z_over_chi[k]*(vars.h[i][k]*d1[j].chi + vars.h[j][k]*d1[i].chi - vars.h[i][j]*d1[k].chi + d1[k].h[i][j]*vars.chi);
         }

         LL[i][j] = (ricci_chi + vars.chi*ricci_tilde + z_terms) / vars.chi;
      }
      scalar = vars.chi*trace(LL, h_UU);
   }

   ricciZ_t(const vars_t<data_t> &vars, const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
            const vars_t<data_t> (&d1)[CH_SPACEDIM],
            const christ_t& chris,
            const tensor<2, data_t>& h_UU)
   {
      tensor<1,double> Z0 = {0};
      ricciZ_t(vars,d2,d1,chris,h_UU,Z0);
   }
};

template <class data_t>
struct ricci_t : ricciZ_t  //Make Z-less ricci explicit for better readibility in the main code.
{
   ricci_t(const vars_t<data_t> &vars, const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
           const vars_t<data_t> (&d1)[CH_SPACEDIM],
           const christ_t& chris,
           const tensor<2, data_t>& h_UU) : ricciZ_t(vars, d2, d1, chris, h_UU) {}
};

template <class data_t>
ALWAYS_INLINE
data_t
compute_trace(const tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &inverse_metric)
{
   data_t trace = 0;
   FOR2(i,j)
   {
      trace += inverse_metric[i][j]*tensor_LL[i][j];
   }
   return trace;
}

template <class data_t>
ALWAYS_INLINE
void
make_trace_free(tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &metric, const tensor<2,data_t> &inverse_metric)
{
   auto trace = compute_trace(tensor_LL, inverse_metric);
   FOR2(i,j)
   {
      tensor_LL[i][j] += - 1./(GR_SPACEDIM-1.) * metric[i][j] * trace;
   }
}

template <class data_t>
ALWAYS_INLINE
tensor<2,data_t>
raise(const tensor<2,data_t> &tensor_LL, const tensor<2,data_t> &inverse_metric)
{
   tensor<2, data_t> tensor_UU;
   FOR2(i,j)
   {
      tensor_UU[i][j] = 0;
      FOR2(k,l)
      {
         tensor_UU[i][j] += inverse_metric[i][k]*inverse_metric[j][l]*tensor_LL[k][l];
      }
   }
   return tensor_UU;
}


