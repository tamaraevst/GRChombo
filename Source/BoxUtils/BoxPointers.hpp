#ifndef BOXPOINTERS_HPP_
#define BOXPOINTERS_HPP_

#include "FArrayBox.H"
#include "UserVariables.hpp"
#include "CellIndex.hpp"

class BoxPointers
{
public:
    const double *m_in_ptr[c_NUM];
    const int *m_in_lo;
    const int *m_in_hi;
    int m_in_stride[3];

    double *m_out_ptr[c_NUM];
    const int *m_out_lo;
    const int *m_out_hi;
    int m_out_stride[3];

    BoxPointers(const FArrayBox& in, FArrayBox& out)
    {
        // dataPtr in Chombo does CH_assert bound check
        // which we don't want to do in a loop
        for (int i = 0; i < c_NUM; ++i) m_in_ptr[i] = in.dataPtr(i);
        //If the output FArrayBox doesn't have all components, just don't set the pointers
        for (int i = 0; i < out.nComp(); ++i) m_out_ptr[i] = out.dataPtr(i);

        m_in_lo = in.loVect();
        m_in_hi = in.hiVect();
        m_in_stride[0] = 1;
        m_in_stride[1] = m_in_hi[0]-m_in_lo[0]+1;
#if CH_SPACEDIM >= 3
        m_in_stride[2] = (m_in_hi[1]-m_in_lo[1]+1)*m_in_stride[1];
#endif
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif

        m_out_lo = out.loVect();
        m_out_hi = out.hiVect();

        m_out_stride[0] = 1;
        m_out_stride[1] = m_out_hi[0]-m_out_lo[0]+1;
#if CH_SPACEDIM >= 3
        m_out_stride[2] = (m_out_hi[1]-m_out_lo[1]+1)*m_out_stride[1];
#endif
#if CH_SPACEDIM >= 4
#error "TODO: Implement CH_SPACEDIM >= 4"
#endif
    }

    CellIndexIn get_in_index(IntVect integer_coords) const
    {
        return ( m_in_stride[2]*(integer_coords[2]-m_in_lo[2]) +
                 m_in_stride[1]*(integer_coords[1]-m_in_lo[1]) + (integer_coords[0]-m_in_lo[0]) );
    }

    CellIndexOut get_out_index(IntVect integer_coords) const
    {
        return ( m_out_stride[2]*(integer_coords[2]-m_out_lo[2]) +
                 m_out_stride[1]*(integer_coords[1]-m_out_lo[1]) + (integer_coords[0]-m_out_lo[0]) );
    }
};

#endif /* BOXPOINTERS_HPP_ */
