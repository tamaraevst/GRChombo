#ifndef GRLEVELDATA_HPP_
#define GRLEVELDATA_HPP_

#include "LevelData.H"
#include "FArrayBox.H"

class GRLevelData : public LevelData<FArrayBox>
{
public:
    void setVal(const double a_val)
    {
        DataIterator dit = m_disjointBoxLayout.dataIterator();
        for(dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& fab = (*this)[dit];
            fab.setVal(a_val); //TODO: Make sure these functions are threaded in Chombo!
        }
    }

    void setVal(const double a_val, const int a_comp)
    {
        DataIterator dit = m_disjointBoxLayout.dataIterator();
        for(dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& fab = (*this)[dit];
            fab.setVal(a_val, a_comp);
        }
    }

    void setVal(const double a_val, const Interval a_comps)
    {
        DataIterator dit = m_disjointBoxLayout.dataIterator();
        //Want component loop inside so unfortunately we have to duplicate the outer loop
        for(dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& fab = (*this)[dit];
            for (int i=a_comps.begin(); i<=a_comps.end(); ++i)
            {
                fab.setVal(a_val, i);
            }
        }
    }

    //a_src and this must have the same box layout
    void plus(const GRLevelData& a_src, const double a_scale)
    {
        DataIterator dit = m_disjointBoxLayout.dataIterator();
        for(dit.begin(); dit.ok(); ++dit)
        {
            FArrayBox& fab           = (*this)[dit];
            const FArrayBox& src_fab = (*this)[dit];
            fab.plus(src_fab, a_scale);
        }
    }
};

#endif /* GRLEVELDATA_HPP_ */
