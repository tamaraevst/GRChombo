#ifndef GRPARMPARSE_HPP_
#define GRPARMPARSE_HPP_

#include "ParmParse.H"
#include <memory>

// Begin: Helper structs to translate a dataype into a Chombo ParmParse data
// type
template <class T> struct ParmParseTranslator;

template <> struct ParmParseTranslator<double>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppDouble;
};

template <> struct ParmParseTranslator<int>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppInt;
};

template <> struct ParmParseTranslator<bool>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppBool;
};
// End: Helper structs to translate a dataype into a Chombo ParmParse data type

class GRParmParse : public ParmParse
{
  public:
    using ParmParse::ParmParse; // Just use ParmParse's constructor

    // Note (MK): I called the functions below "load" rather than "get" to avoid
    // clashes with the many  different overloads of "get" in ParmParse. Also, I
    // think load is a more intuitive name.

    /// Loads an array from the parameter file
    template <class data_t, long unsigned int n_comp>
    void load(const char *name, std::array<data_t, n_comp> &array) const
    {
        std::unique_ptr<data_t[]> c_array{new data_t[n_comp]};
        getarr(name, ParmParseTranslator<data_t>::pp_type, c_array.get(), 0,
               n_comp, -1);
        for (int i = 0; i < n_comp; ++i)
            array[i] = c_array[i];
    }

    /// Loads a vector with num_comp components from the parameter file
    template <class data_t>
    void load(const char *name, std::vector<data_t> &vector,
              const int num_comp) const
    {
        getarr(name, vector, 0, num_comp);
    }

    /// Loads a value from the paramter file
    template <class data_t>
    typename std::enable_if<
        !std::is_enum<data_t>::value>::type // Can't use for enum types
    load(const char *name, data_t &parameter) const
    {
        get(name, parameter);
    }

    /// Loads an enum value from the paramter file
    template <typename enum_type>
    typename std::enable_if<
        std::is_enum<enum_type>::value>::type // Only enabled for enum types
    load(const char *name, enum_type &parameter) const
    {
        int iparam;
        get(name, iparam);
        parameter = static_cast<enum_type>(iparam);
    }

    /// Loads a value from the paramter file, if the value isn't defined it sets
    /// to the supplied default
    template <class data_t>
    void load(const char *name, data_t &parameter,
              const data_t default_value) const
    {
        if (contains(name))
        {
            load(name, parameter);
        }
        else
        {
            parameter = default_value;
            pout() << "Parameter: " << name << "not found in parameter file. "
                   << "It has been set to its default value." << std::endl;
        }
    }
};

#endif /* GRPARMPARSE_HPP_ */
