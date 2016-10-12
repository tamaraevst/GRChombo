#ifndef SFMATTER_HPP_
#define SFMATTER_HPP_

#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.hpp"
#include "FABDriverBase.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "CCZ4Geometry.hpp"
#include "VarsBase.hpp"
#include "CCZ4.hpp"

#include "UserVariables.hpp" //This files needs c_NUM - total number of components

#include <array>

class SFMatter
{

public:

    SFMatter(void)
    {};

    template <class data_t>
    struct vars_t : VarsBase<data_t>
    {
        using VarsBase<data_t>::define_enum_mapping; //Saves us some writing later
        using VarsBase<data_t>::define_symmetric_enum_mapping;

        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;
				data_t phi;
        data_t Pi;

        vars_t();
    };

//    template <class data_t>
//    struct matter_vars_t
//    {
//				data_t matter_phi;
//        data_t matter_Pi;
//    };

    template <class data_t>
		struct emtensor_t
		{
      tensor<2, data_t> Sij; // S_ij = T_ij
      tensor<1,data_t>  Si; // S_i = T_ia_n^a
      data_t            S; // S = S^i_i
      data_t            rho; // rho = T_ab n^a n^b
		};

    template <class data_t>
		struct potential_t
		{
      data_t						V_of_phi; //V(\phi)
      data_t            dVdphi; // Gradient of V(\phi)
		};

    template <class data_t>
    emtensor_t<data_t> calc_emtensor(
         const vars_t<data_t> &vars,
         const vars_t< tensor<1,data_t> >& d1,
         const tensor<2, data_t>& h_UU,
         const tensor<3, data_t>& chris,
         const vars_t<data_t> &advec
    );

    template <class data_t>
    potential_t<data_t> calc_potential(
				 const data_t phi
    );

		template <class data_t>
		vars_t<data_t> calc_total_rhs(
        const CCZ4::vars_t<data_t> &CCZ4_rhs,
        const vars_t<data_t> &matter_rhs,
        const vars_t<data_t> &vars,
        const vars_t< tensor<1,data_t> >& d1,
        const vars_t< tensor<2,data_t> >& d2,
        const vars_t<data_t> &advec
    );

};

#include "SFMatter.impl.hpp"

#endif /* SFMATTER_HPP_ */
