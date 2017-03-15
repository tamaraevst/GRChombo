#ifndef _ARRAY_HPP_
#define _ARRAY_HPP_

// Basically implementing the basic capabilities of std::array
template <typename T, int N>
class Array
{
protected:
    T arr[N];

    inline void
    copyFrom(const T (&rhs_arr)[N])
    {
        for (int i = 0; i < N; ++i)
        {
            arr[i] = rhs_arr[i];
        }
    }

public:
    inline
    Array()
    {

    }

    inline
    Array(const Array<T, N>& other)
    {
        copyFrom(other.arr);
    }

    inline
    Array(const T (&other)[N])
    {
        copyFrom(other);
    }

    inline Array<T, N>&
    operator=(const Array<T, N>& rhs)
    {
        copyFrom(rhs.arr);
        return *this;
    }

    inline Array<T, N>&
    operator=(const T (&other)[N])
    {
        copyFrom(other);
        return *this;
    }

    inline T&
    operator[](int i)
    {
        return arr[i];
    }

    inline const T&
    operator[](int i) const
    {
        return arr[i];
    }

};



#endif /* _ARRAY_HPP_ */
