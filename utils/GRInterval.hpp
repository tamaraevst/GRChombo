#pragma once

///A templated version of Chombo's Interval - allows compile time checking.
template<int ibegin, int iend>
struct GRInterval
{
    static constexpr int begin()
    {
        return ibegin;
    }

    static constexpr int end()
    {
        return iend;
    }

    static constexpr int size()
    {
        return iend-ibegin;
    }

    static constexpr bool contains(int i)
    {
        return ((i>ibegin) && (i<iend));
    }
};
