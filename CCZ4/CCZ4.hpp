#ifndef CCZ4_HPP_
#define CCZ4_HPP_

#include "FArrayBox.H"
#include "simd.hpp"
#include "tensor.hpp"
#include "GRUtils.H"
#include "FABDriverBase.hpp"

#include "user_enum.hpp" //TODO: Check whether this is really necessary!

class CCZ4
{
public:
    template <class data_t>
    struct vars_t
    {
        data_t chi;
        tensor<2, data_t> h;
        data_t K;
        tensor<2, data_t> A;
        tensor<1, data_t> Gamma;
        data_t Theta;
        data_t lapse;
        tensor<1, data_t> shift;
        tensor<1, data_t> B;
    };

public:
    struct params_t
    {
        double kappa1;
        double kappa2;
        double kappa3;
        double shift_gamma_coeff;
        double lapse_advec_coeff;
        double shift_advec_coeff;
        double beta_driver;
    };

protected:
    params_t m_params;
    double m_dx;
    double m_sigma;

//    const double *m_in_ptr[c_NUM];  //TODO: Delete all commented out bits once it compiles
//    const int *m_in_lo;
//    const int *m_in_hi;
//    int m_stride[3];
//
//    double *m_out_ptr[c_NUM];
//    const int *m_out_lo;
//    const int *m_out_hi;
//    int m_out_stride[3];
//
public:
    CCZ4(params_t params, double dx, double sigma, const FABDriverBase& driver);
    void execute(const FArrayBox& in, FArrayBox& out);

//    class factory
//    {
//    public:
//        params_t m_params;
//        double m_dx;
//        double m_sigma;
//
//        factory(params_t params, double dx, double sigma) :
//            m_params (params),
//            m_dx (dx),
//            m_sigma (sigma)
//        {}
//
//        CCZ4 create(FABDriver &driver)
//        {
//            return CCZ4(params, dx, sigma, driver);
//        }
//    }

protected:
    const FABDriverBase& m_driver;

    template <class data_t>
    void compute(int x, int y, int z);

    template <class data_t>
    void demarshall(const data_t (&in)[c_NUM], vars_t<data_t>& out);

    template <class data_t>
    void rhs_equation(const vars_t<data_t> &vars,
             const vars_t<data_t> (&d1)[CH_SPACEDIM],
             const vars_t<data_t> (&d2)[CH_SPACEDIM][CH_SPACEDIM],
             const vars_t<data_t> &advec,
             vars_t<data_t> &rhs);
};
#endif /* CCZ4_HPP_ */
