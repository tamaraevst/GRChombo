#ifndef NANCHECK_HPP_
#define NANCHECK_HPP_

// This class offers a nan-check including debugging information when it happens

#include "Cell.hpp"
#include "UserVariables.hpp"

class NanCheck
{
  protected:
    const std::string m_error_info = "NanCheck";
    const double m_max_abs = 1e20;

  public:
    NanCheck() {}

    NanCheck(const std::string a_error_info) : m_error_info(a_error_info) {}

    NanCheck(const std::string a_error_info, const double a_max_abs)
        : m_error_info(a_error_info), m_max_abs(a_max_abs)
    {
    }

    void compute(Cell<double> current_cell) const
    {
        bool stop = false;
        for (int ivar = 0; ivar < c_NUM; ++ivar)
        {
            double val;
            current_cell.load_vars(val, ivar);
            if (std::isnan(val) || abs(val) > m_max_abs)
                stop = true;
        }
        if (stop)
        {
#pragma omp single
            {
                pout() << m_error_info
                       << "::Values have become nan. The current state is: "
                       << endl;
                for (int ivar = 0; ivar < c_NUM; ++ivar)
                {
                    pout() << UserVariables::variable_names[ivar] << ": "
                        << current_cell.load_vars(ivar) << endl;
                }
                pout() << "Integer coordinates: " << current_cell.get_int_vect()
                       << endl;
            }
            MayDay::Error("Values have become nan.");
        }
    }
};

#endif /* NANCHECK_HPP_ */
