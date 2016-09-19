#ifndef STACKVECTOR_HPP_
#define STACKVECTOR_HPP_

template <class T, int N>
class StackVector
{
private:
    int  m_ncomp = 0;

public:
    T    m_data[N];

    int
    push_back(T new_datum)
    {
        if (m_ncomp < N)
        {
            m_data[m_ncomp] = new_datum;
            m_ncomp++;
            return 0;
        }
        else
            return -1;
    }

    T& operator [](int i)
    {
        return m_data[i];
    }

    const T& operator [](int i) const
    {
        return m_data[i];
    }

    int
    get_ncomp() const
    {
        return m_ncomp;
    }
};
#endif /* STACKVECTOR_HPP_ */
