#ifndef INDEX_APPLICATOR_HPP
#define INDEX_APPLICATOR_HPP

//This class allows one to apply an arbitrary number of indices to an object
//of an arbitrary number of indices.
//This is useful for example in VarsBase where the IndexApplicator allows us to write general
//code for an arbitrary number of derivative indices.

//Need a base template that we can specialise. Note that this will never be instantiated
//because the specialisations cover everything.
template <typename... Ts>
class IndexApplicator
{
};

//Specialisation for one or more indices
template <typename t, typename... Ts>
class IndexApplicator<t, Ts...>
{
public:
    template <typename data_t>
    static ALWAYS_INLINE auto //Always inlined so this iteration trick doesn't ruin performance
    apply(data_t& obj, t dir0, Ts... dirs)
    -> decltype(IndexApplicator<Ts...>::apply(obj[dir0], dirs...))&
    {
        //Let the compiler iterate until there are no indices left (at which point we go to the
        //specialisation below).
        return IndexApplicator<Ts...>::apply(obj[dir0], dirs...);
    }
};

//Specialisation for the case with no indices left
template <>
class IndexApplicator<>
{
public:
    template <typename data_t>
    static ALWAYS_INLINE data_t&
    apply(data_t& obj)
    {
        return obj;
    }
};

#endif /* INDEX_APPLICATOR_HPP */
