#include <iostream>
#include "FABDriverBase.hpp"
#include "VarsBase.hpp"

struct vars_t : VarsBase<double>
{
    double var;
    double symmetric_var_1;
    double symmetric_var_2;

    vars_t()
    {
        define_enum_mapping(c_var, var);
        define_enum_mapping(c_sym_var, symmetric_var_1);
        define_enum_mapping(c_sym_var, symmetric_var_2);
    }
};

int main()
{
    int failed = 0;

    vars_t vars;
    vars.var = 42.;
    vars.symmetric_var_1 = 84.;
    vars.symmetric_var_2 = 84.;


    double out[c_NUM];
    FABDriverBase driver;
    FORVARS(i) driver.m_out_ptr[i] = &out[i];

    driver.store_vars(vars);

    if (out[c_var] != 42.) failed=1;
    if (out[c_sym_var] != 84.) failed=2;

    if (failed) std::cout << "Variable store test FAILED with code " << failed << std::endl;
    else std::cout << "Variable store test PASSED. (return code " << failed << ")" << std::endl;
    return failed;
}
