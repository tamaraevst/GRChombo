#include <iostream>

#include "UserVariables.hpp"
#include "Cell.hpp"
#include "VarsBase.hpp"
#include "CellIndex.hpp"

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


    Box box(IntVect(0,0,0), IntVect(0,0,0));
    FArrayBox fab_in(box, 3);
    FArrayBox fab_out(box, 3);
    Cell<data_t> current_cell(IntVect(0,0,0), BoxPointers(fab_in, fab_out));

    current_cell.store_vars(vars);

    if (fab_out(IntVect::Zero,0) != 42.) failed=1;
    if (fab_out(IntVect::Zero,1) != 84.) failed=2;

    if (failed) std::cout << "Variable store test FAILED with code " << failed << std::endl;
    else std::cout << "Variable store test PASSED. (return code " << failed << ")" << std::endl;
    return failed;
}
