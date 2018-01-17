#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

enum
{
    c_var,
    c_sym_var,
    c_NUM
};

namespace UserVariables
{
static constexpr char const *variable_names[c_NUM] = {
    "var", "sym_var",
};
}

#endif /* USERVARIABLES_HPP_ */
