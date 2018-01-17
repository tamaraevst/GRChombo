#pragma once

#include "GRInterval.hpp"

namespace VarsTools
{
template <typename mapping_function_t, typename data_t>
void define_enum_mapping(mapping_function_t mapping_function, const int &ivar,
                         data_t &scalar)
{
    mapping_function(ivar, scalar);
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
void define_enum_mapping(mapping_function_t mapping_function,
                         const GRInterval<start_var, end_var> interval,
                         tensor<1, data_t, end_var - start_var + 1> &tensor)
{
    for (int ivar = 0; ivar < interval.size(); ++ivar)
        mapping_function(start_var + ivar, tensor[ivar]);
}

template <typename mapping_function_t, typename data_t, int start_var,
          int end_var>
void define_symmetric_enum_mapping(
    mapping_function_t mapping_function,
    const GRInterval<start_var, end_var> interval, tensor<2, data_t> &tensor)
{
    static_assert(interval.size() == IDX_SPACEDIM * (IDX_SPACEDIM + 1) / 2,
                  "Interval has wrong size");
#if IDX_SPACEDIM == 3
    mapping_function(start_var, tensor[0][0]);

    mapping_function(start_var + 1, tensor[0][1]);
    mapping_function(start_var + 1, tensor[1][0]);

    mapping_function(start_var + 2, tensor[0][2]);
    mapping_function(start_var + 2, tensor[2][0]);

    mapping_function(start_var + 3, tensor[1][1]);

    mapping_function(start_var + 4, tensor[1][2]);
    mapping_function(start_var + 4, tensor[2][1]);

    mapping_function(start_var + 5, tensor[2][2]);
#else
#error IDX_SPACEDIM not equal to three not implemented yet...
#endif
}

//--> Begin: Helper for the assign function
template <class nested_template> struct _strip_nested_template;

template <template <typename> class outermost_layer, class inner_part>
struct _strip_nested_template<outermost_layer<inner_part>>
{
    using type = inner_part;
};
//<-- End: Helper for the assign function

/// Writes data directly into all variables
/**if this variables has multiple components (e.g. if it is an array of
 *derivatives) the data can be written directly into these components by
 *specifying an arbitrary number of icomps
 */
template <class vars_t, typename value_t>
ALWAYS_INLINE void assign(vars_t &vars, const value_t &value)
{
    // The template magic below is needed to make sure that we can write
    // assign(vars, 0.)  and 0. gets correctly cast from double to simd<double> if
    // necessary.
    using data_t = typename _strip_nested_template<vars_t>::type;
    vars.enum_mapping([&value](const int &ivar, data_t &var) {
        var = static_cast<data_t>(value);
    });
}

/// Prints all elements of the vars element with component names - Useful for
/// debugging.
template <template <typename> class vars_t, typename data_t>
void print(const vars_t<data_t> &vars)
{
    vars.enum_mapping([](const int &ivar, data_t &var) {
        pout() << UserVariables::variable_names[ivar] << ": " << var << "\n";
    });
}
}
