#ifndef USERVARIABLES_HPP_
#define USERVARIABLES_HPP_

enum {
    c_chi,
    c_lapse,
    c_NUM
};

namespace UserVariables
{
    static constexpr char const * variable_names[c_NUM] =
    {
        "chi",
        "lapse",
    };
}

#endif /* USERVARIABLES_HPP_ */
